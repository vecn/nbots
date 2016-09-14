#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "nb/math_bot.h"
#include "nb/memory_bot.h"
#include "nb/interpolation_bot.h"
#include "nb/eigen_bot.h"
#include "nb/geometric_bot.h"
#include "nb/graph_bot.h"
#include "nb/pde_bot/material.h"
#include "nb/pde_bot/gauss_legendre_quad.h"
#include "nb/pde_bot/common_solid_mechanics/analysis2D.h"
#include "nb/pde_bot/common_solid_mechanics/formulas.h"
#include "nb/pde_bot/boundary_conditions/bcond.h"
#include "nb/pde_bot/boundary_conditions/bcond_iter.h"

#include "set_bconditions.h"

#define QUADRATURE_POINTS 4
#define ALL_NEIGHBOURS true


#define POW2(a) ((a)*(a))

static int assemble_system(vcn_sparse_t *K, double *F,
			   const nb_partition_t *const part,
			   const nb_material_t *material,
			   bool enable_self_weight,
			   double gravity[2],
			   nb_analysis2D_t analysis2D,
			   nb_analysis2D_params *params2D);
static void assemble_elem(uint32_t elem_id,
			  vcn_sparse_t *K, double *F,
			  const nb_partition_t *const part,
			  const nb_material_t *material,
			  bool enable_self_weight,
			  double gravity[2],
			  nb_analysis2D_t analysis2D,
			  nb_analysis2D_params *params2D);
static void integrate_elem_force(const nb_partition_t *part,
				 const nb_material_t *material,
				 bool enable_self_weight,
				 double gravity[2],
				 uint32_t elem_id,
				 double *F);
static void assemble_face(uint32_t elem_id, uint16_t face_id,
			  const nb_partition_t *const part,
			  const nb_material_t *material,
			  nb_analysis2D_t analysis2D,
			  nb_analysis2D_params *params2D,
			  vcn_sparse_t *K,  double *F);
static void integrate_face(uint32_t elem_id, uint16_t face_id,
			   const nb_partition_t *const part,
			   const nb_material_t *material,
			   nb_analysis2D_t analysis2D,
			   nb_analysis2D_params *params2D,
			   vcn_sparse_t *K);
static uint16_t get_neighbours(const nb_partition_t *part, uint32_t elem_id,
			       uint16_t face_id, uint32_t *ngb);
static uint16_t get_ngb_around_right_vtx(const nb_partition_t *part,
					 uint32_t *ngb, uint16_t current_id,
					 uint32_t elem_id,
					 uint16_t face_id);
static void get_Kf(const nb_partition_t *const part,
		   uint32_t elem_id, uint16_t face_id,
		   const double D[4], uint16_t N_ngb,
		   const uint32_t *ngb, uint8_t N_qp,
		   double *Kf);
static void get_quadrature_points(const nb_partition_t *part,
				  uint32_t elem_id, uint16_t face_id,
				  double lf, uint8_t N_qp,
				  double *xqp, double *wqp);
static void interpolators_eval_grad(const nb_partition_t *part, uint8_t N_ngb,
				    const uint32_t *ngb, const double x[2],
				    double *grad_phi);
static void get_Kf_nodal_contribution(const nb_partition_t *part,
				      const double D[4], const double nf[2],
				      uint16_t i, const double *grad_phi,
				      double Kfi[4]);
static void add_Kf_to_K(vcn_sparse_t *K, const double *Kf,
			uint32_t elem_id, uint32_t ngb_id,
			uint16_t N_ngb, const uint32_t *ngb,
			double factor);
static int solver(const vcn_sparse_t *const A,
		  const double *const b, double* x);
static void compute_strain(double *strain,
			   const nb_partition_t *const part,
			   double *disp,
			   nb_analysis2D_t analysis2D,
			   const nb_material_t *const material);

int nb_cvfa_compute_2D_Solid_Mechanics
			(const nb_partition_t *const part,
			 const nb_material_t *const material,
			 const nb_bcond_t *const bcond,
			 bool enable_self_weight,
			 double gravity[2],
			 nb_analysis2D_t analysis2D,
			 nb_analysis2D_params *params2D,
			 double *displacement, /* Output */
			 double *strain       /* Output */)
{
	int status = 0;
	vcn_graph_t *graph = malloc(nb_graph_get_memsize());
	nb_graph_init(graph);
	nb_partition_load_graph(part, graph, NB_ELEMS_LINKED_BY_NODES);
	vcn_sparse_t *K = vcn_sparse_create(graph, NULL, 2);
	nb_graph_finish(graph);

	uint32_t N_elems = nb_partition_get_N_elems(part);
	uint32_t F_memsize = 2 * N_elems * sizeof(double);
	double* F = NB_SOFT_MALLOC(F_memsize);
	memset(F, 0, F_memsize);

	int status_assemble = assemble_system(K, F, part, material,
					      enable_self_weight, gravity,
					      analysis2D, params2D);


	if (0 != status_assemble) {
		status = 1;
		goto CLEANUP_LINEAR_SYSTEM;
	}

	nb_cvfa_set_bconditions(part, K, F, bcond, 1.0);

	int solver_status = solver(K, F, displacement);
	if (0 != solver_status) {
		status = 2;
		goto CLEANUP_LINEAR_SYSTEM;
	}

	compute_strain(strain, part, displacement,
		       analysis2D, material);

	status = 0;
CLEANUP_LINEAR_SYSTEM:
	vcn_sparse_destroy(K);
	NB_SOFT_FREE(F_memsize, F);
	return status;
}

static int assemble_system(vcn_sparse_t *K, double *F,
			   const nb_partition_t *const part,
			   const nb_material_t *material,
			   bool enable_self_weight,
			   double gravity[2],
			   nb_analysis2D_t analysis2D,
			   nb_analysis2D_params *params2D)
{
	uint32_t N_elems = nb_partition_get_N_elems(part);
	vcn_sparse_reset(K);
	memset(F, 0, vcn_sparse_get_size(K) * sizeof(*F));

	for (uint32_t i = 0; i < N_elems; i++) {
		assemble_elem(i, K, F, part, material,
			      enable_self_weight, gravity,
			      analysis2D, params2D);
			
	}

	return 0;
}

static void assemble_elem(uint32_t elem_id,
			  vcn_sparse_t *K, double *F,
			  const nb_partition_t *const part,
			  const nb_material_t *material,
			  bool enable_self_weight,
			  double gravity[2],
			  nb_analysis2D_t analysis2D,
			  nb_analysis2D_params *params2D)
{
	integrate_elem_force(part, material, enable_self_weight,
			     gravity, elem_id, F);

	uint16_t N_faces = nb_partition_elem_get_N_adj(part, elem_id);
	for (uint16_t j = 0; j < N_faces; j++) {
		assemble_face(elem_id, j, part, material,
			      analysis2D, params2D,
			      K, F);
	}
}

static void integrate_elem_force(const nb_partition_t *part,
				 const nb_material_t *material,
				 bool enable_self_weight,
				 double gravity[2],
				 uint32_t elem_id,
				 double *F)
{
	if (enable_self_weight) {
		double area = nb_partition_elem_get_area(part, elem_id);
		double mass = area * nb_material_get_density(material);
		F[elem_id * 2] += mass * gravity[0];
		F[elem_id*2+1] += mass * gravity[1];
	}
}

static void assemble_face(uint32_t elem_id, uint16_t face_id,
			  const nb_partition_t *const part,
			  const nb_material_t *material,
			  nb_analysis2D_t analysis2D,
			  nb_analysis2D_params *params2D,
			  vcn_sparse_t *K, double *F)
{
	if (nb_partition_elem_has_ngb(part, elem_id, face_id)) {
		uint32_t ngb_id = nb_partition_elem_get_ngb(part, elem_id,
							    face_id);
		if (elem_id < ngb_id)
			/* Integrate face of both elements */
			integrate_face(elem_id, face_id, part, material,
				       analysis2D, params2D, K);
	}
}

static void integrate_face(uint32_t elem_id, uint16_t face_id,
			   const nb_partition_t *const part,
			   const nb_material_t *material,
			   nb_analysis2D_t analysis2D,
			   nb_analysis2D_params *params2D,
			   vcn_sparse_t *K)
{
	double D[4];
	nb_pde_get_constitutive_matrix(D, material, analysis2D);

	uint32_t ngb[16];
	uint16_t N_ngb = get_neighbours(part, elem_id, face_id, ngb);
	
	uint32_t memsize = 4 * N_ngb * sizeof(double);
	double *Kf = NB_SOFT_MALLOC(memsize);
	get_Kf(part, elem_id, face_id, D, N_ngb, ngb, QUADRATURE_POINTS, Kf);

	double factor = params2D->thickness;

	uint32_t ngb_id = nb_partition_elem_get_ngb(part, elem_id, face_id);
	add_Kf_to_K(K, Kf, elem_id, ngb_id, N_ngb, ngb, factor);

	NB_SOFT_FREE(memsize, Kf);
}

static uint16_t get_neighbours(const nb_partition_t *part, uint32_t elem_id,
			       uint16_t face_id, uint32_t *ngb)
{
	ngb[0] = elem_id;
	ngb[1] = nb_partition_elem_get_ngb(part, elem_id, face_id);
	uint16_t N = 2;

	if (ALL_NEIGHBOURS) {
		N = get_ngb_around_right_vtx(part, ngb, N, ngb[0], face_id);
		uint16_t aux = nb_partition_elem_ngb_get_face(part, ngb[1],
							      ngb[0]);
		N = get_ngb_around_right_vtx(part, ngb, N, ngb[1], aux);
	}
	return N;
}

static uint16_t get_ngb_around_right_vtx(const nb_partition_t *part,
					 uint32_t *ngb, uint16_t current_id,
					 uint32_t elem_id,
					 uint16_t face_id)
{
	uint32_t N_elems = nb_partition_get_N_elems(part);

	uint32_t front_ngb_id = nb_partition_elem_get_ngb(part, elem_id, face_id);
	uint32_t nid_prev = elem_id;
	uint32_t nid = nb_partition_elem_face_get_right_ngb(part, elem_id, face_id);
	while (nid != front_ngb_id && nid < N_elems) {
		ngb[current_id] = nid;
		current_id += 1;
		uint16_t aux = nb_partition_elem_ngb_get_face(part, nid,
							      nid_prev);
		nid_prev = nid;
		nid = nb_partition_elem_face_get_right_ngb(part, nid, aux);
	}
	if (nid >= N_elems && front_ngb_id < N_elems) {
		nid_prev = front_ngb_id;
		uint16_t aux = nb_partition_elem_ngb_get_face(part, front_ngb_id,
							      elem_id);
		nid = nb_partition_elem_face_get_left_ngb(part, front_ngb_id, aux);
		while (nid < N_elems) {
			ngb[current_id] = nid;
			current_id += 1;
			aux = nb_partition_elem_ngb_get_face(part, nid,
							     nid_prev);
			nid_prev = nid;
			nid = nb_partition_elem_face_get_left_ngb(part, nid,
								  aux);
		}
	}
	return current_id;
}

static void get_Kf(const nb_partition_t *const part,
		   uint32_t elem_id, uint16_t face_id,
		   const double D[4], uint16_t N_ngb,
		   const uint32_t *ngb, uint8_t N_qp,
		   double *Kf)
{
	uint32_t memsize = (3 * N_qp + 2 * N_ngb) * sizeof(double);
	char *memblock = NB_SOFT_MALLOC(memsize);
	double *wqp = (void*) memblock;
	double *xqp = (void*) (memblock + N_qp * sizeof(double));
	double *grad_phi = (void*) (memblock + 3 * N_qp * sizeof(double));

	double nf[2];
	double lf = nb_partition_elem_face_get_normal(part, elem_id,
						      face_id, nf);

	get_quadrature_points(part, elem_id, face_id, lf, N_qp, xqp, wqp);

	memset(Kf, 0, 4 * N_ngb * sizeof(double));
	for (uint8_t q = 0; q < N_qp; q++) {
		interpolators_eval_grad(part, N_ngb, ngb, 
					&(xqp[q*2]), grad_phi);
	
		for (uint16_t i = 0; i < N_ngb; i++) {
			double Kfi[4];
			get_Kf_nodal_contribution(part, D, nf, i,
						  grad_phi, Kfi);
			Kf[i * 2] += wqp[q] * Kfi[0];
			Kf[i*2+1] += wqp[q] * Kfi[1];
			Kf[2 * N_ngb + i * 2] += wqp[q] * Kfi[2];
			Kf[2 * N_ngb + i*2+1] += wqp[q] * Kfi[3];
		}
	}
	NB_SOFT_FREE(memsize, memblock);
}

static void get_quadrature_points(const nb_partition_t *part,
				  uint32_t elem_id, uint16_t face_id,
				  double lf, uint8_t N_qp,
				  double *xqp, double *wqp)
{
	nb_glquadrature_t glq;
	glq.x = alloca(N_qp * sizeof(*(glq.x)));
	glq.w = alloca(N_qp * sizeof(*(glq.w)));
	nb_glquadrature_load(&glq, N_qp);

	for (uint8_t q = 0; q < N_qp; q++) {
		wqp[q] = 0.5 * glq.w[q] * lf;
		double w = (glq.x[q] + 1)/2.0;
		nb_partition_elem_face_get_midpoint(part, elem_id, face_id,
						    w, &(xqp[q*2]));
	}
}

static void interpolators_eval_grad(const nb_partition_t *part, uint8_t N_ngb,
				    const uint32_t *ngb, const double x[2],
				    double *grad_phi)
{	
	uint32_t memsize = 2 * N_ngb * sizeof(double);
	double *ni = NB_SOFT_MALLOC(memsize);

	for (uint32_t i = 0; i < N_ngb; i++) {
		ni[i * 2] = nb_partition_elem_get_x(part, ngb[i]);
		ni[i*2+1] = nb_partition_elem_get_y(part, ngb[i]);
	}

	nb_nonpolynomial_simple_eval_grad(N_ngb, 2, ni, x, grad_phi);	

	NB_SOFT_FREE(memsize, ni);
}

static void get_Kf_nodal_contribution(const nb_partition_t *part,
				      const double D[4], const double nf[2],
				      uint16_t i, const double *grad_phi,
				      double Kfi[4])
{
	double dphi_dx = grad_phi[i * 2];
	double dphi_dy = grad_phi[i*2+1];
	Kfi[0] = nf[0] * D[0] * dphi_dx + nf[1] * D[3] * dphi_dy;
	Kfi[1] = nf[0] * D[1] * dphi_dy + nf[1] * D[3] * dphi_dx;
	Kfi[2] = nf[1] * D[1] * dphi_dx + nf[0] * D[3] * dphi_dy;
	Kfi[3] = nf[1] * D[2] * dphi_dy + nf[0] * D[3] * dphi_dx;
}

static void add_Kf_to_K(vcn_sparse_t *K, const double *Kf,
			uint32_t elem_id, uint32_t ngb_id,
			uint16_t N_ngb, const uint32_t *ngb,
			double factor)
{
	uint16_t size = 2 * N_ngb;
	factor *= -1;

	/* Add equations of element */
	uint32_t i = elem_id;
	for (uint32_t k = 0; k < N_ngb; k++) {
		uint32_t j = ngb[k];
		uint8_t c1 = k * 2;
		uint8_t c2 = k * 2 + 1;
		vcn_sparse_add(K, i * 2, j * 2, factor * Kf[c1]);
		vcn_sparse_add(K, i * 2, j*2+1, factor * Kf[c2]);
		vcn_sparse_add(K, i*2+1, j * 2, factor * Kf[size + c1]);
		vcn_sparse_add(K, i*2+1, j*2+1, factor * Kf[size + c2]);
	}
	/* Add equations of neighbour */
	i = ngb_id;
	for (uint32_t k = 0; k < N_ngb; k++) {
		uint32_t j = ngb[k];
		uint8_t c1 = k * 2;
		uint8_t c2 = k * 2 + 1;
		vcn_sparse_add(K, i * 2, j * 2, -factor * Kf[c1]);
		vcn_sparse_add(K, i * 2, j*2+1, -factor * Kf[c2]);
		vcn_sparse_add(K, i*2+1, j * 2, -factor * Kf[size + c1]);
		vcn_sparse_add(K, i*2+1, j*2+1, -factor * Kf[size + c2]);
	}
}

static int solver(const vcn_sparse_t *const A,
		  const double *const b, double* x)
{
	uint32_t N = vcn_sparse_get_size(A);
	int status = vcn_sparse_solve_using_LU(A, b, x, 1);
	return status;
}

static void compute_strain(double *strain,
			   const nb_partition_t *const part,
			   double *disp,
			   nb_analysis2D_t analysis2D,
			   const nb_material_t *const material)
{
	uint32_t N_elems = nb_partition_get_N_elems(part);
	for (uint32_t i = 0; i < N_elems; i++) {
		double Du[4] = {1,2,3,4};
		strain[i * 3] = Du[0];
		strain[i*3+1] = Du[3];
		strain[i*3+2] = (Du[1] + Du[2]);
	}
}

void nb_cvfa_compute_stress_from_strain(const nb_partition_t *part,
					const nb_material_t *const material,
					nb_analysis2D_t analysis2D,
					double* strain,
					double* stress /* Output */)
{
	uint32_t N_elems = nb_partition_get_N_elems(part);
	for (uint32_t i = 0; i < N_elems; i++) {
		double D[4];
		nb_pde_get_constitutive_matrix(D, material, analysis2D);

		stress[i * 3] = strain[i * 3] * D[0] +
			strain[i*3+1] * D[1];
		stress[i*3+1] = strain[i * 3] * D[1] +
			strain[i*3+1] * D[2];
		stress[i*3+2] = strain[i*3+2] * D[3];
	}
}
