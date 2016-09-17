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

#define QUADRATURE_POINTS 2
#define ALL_NEIGHBOURS true
#define ENABLE_LEAST_SQUARES false


#define POW2(a) ((a)*(a))

static int assemble_system(vcn_sparse_t *K, double *F,
			   const nb_partition_t *const part,
			   const nb_graph_t *graph,
			   const nb_material_t *material,
			   bool enable_self_weight,
			   double gravity[2],
			   nb_analysis2D_t analysis2D,
			   nb_analysis2D_params *params2D);
static void assemble_elem(uint32_t elem_id,
			  vcn_sparse_t *K, double *F,
			  const nb_partition_t *const part,
			  const nb_graph_t *graph,
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
static void integrate_elem(uint32_t elem_id,
			   double *Ke, double *F,
			   const nb_partition_t *const part,
			   const nb_graph_t *graph,
			   const nb_material_t *material,
			   nb_analysis2D_t analysis2D,
			   nb_analysis2D_params *params2D);
static void integrate_face(uint32_t elem_id,
			   uint16_t face_id,
			   const nb_partition_t *const part,
			   const nb_graph_t *graph,
			   const nb_material_t *material,
			   nb_analysis2D_t analysis2D,
			   nb_analysis2D_params *params2D,
			   double *Ke);
static uint16_t face_get_N_ngb(const nb_partition_t *part,
			       uint32_t elem_id, uint16_t face_id);
static uint16_t get_N_ngb_around_right_vtx(const nb_partition_t *part,
					   uint32_t elem_id,
					   uint16_t face_id);
static void face_get_neighbourhood(const nb_partition_t *part,
				   uint32_t elem_id, uint16_t face_id,
				   uint32_t *ngb);
static uint16_t get_ngb_around_right_vtx(const nb_partition_t *part,
					 uint32_t *ngb, uint16_t current_id,
					 uint32_t elem_id,
					 uint16_t face_id);
static void get_Kf(const nb_partition_t *const part,
		   uint32_t elem_id, uint16_t face_id,
		   const double D[4], uint16_t N_ngb,
		   const uint32_t *ngb, uint8_t N_qp,
		   nb_analysis2D_params *params2D,
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
static void add_Kf_to_Ke(double *Ke, const double *Kf,
			 const nb_graph_t *graph, uint32_t elem_id,
			 uint16_t N_ngb, const uint32_t *ngb);
static uint16_t get_id_local_from_global(const nb_graph_t *graph,
					 uint32_t elem_id,
					 uint32_t global_id);
static double get_KTKe(uint16_t N, const double *Ke, double *KTK);
static void add_KTKe_to_K(const nb_graph_t *graph,  uint32_t elem_id,
			  const double *KTK, vcn_sparse_t *K);
static uint32_t get_id_global_from_local(const nb_graph_t *graph,
					 uint32_t elem_id,
					 uint16_t local_id);
static void add_Ke_to_K(const nb_graph_t *graph,  uint32_t elem_id,
			  const double *Ke, vcn_sparse_t *K);
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
	nb_graph_extend_adj(graph, 1);
	vcn_sparse_t *K = vcn_sparse_create(graph, NULL, 2);
	nb_graph_finish(graph);
	nb_graph_init(graph);
	nb_partition_load_graph(part, graph, NB_ELEMS_LINKED_BY_NODES);

	uint32_t N_elems = nb_partition_get_N_elems(part);
	uint32_t F_memsize = 2 * N_elems * sizeof(double);
	double* F = NB_SOFT_MALLOC(F_memsize);
	memset(F, 0, F_memsize);

	int status_assemble = assemble_system(K, F, part, graph, material,
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
	nb_graph_finish(graph);
	vcn_sparse_destroy(K);
	NB_SOFT_FREE(F_memsize, F);
	return status;
}

static int assemble_system(vcn_sparse_t *K, double *F,
			   const nb_partition_t *const part,
			   const nb_graph_t *graph,
			   const nb_material_t *material,
			   bool enable_self_weight,
			   double gravity[2],
			   nb_analysis2D_t analysis2D,
			   nb_analysis2D_params *params2D)
{
	uint32_t N_elems = nb_partition_get_N_elems(part);
	vcn_sparse_reset(K);
	memset(F, 0, vcn_sparse_get_size(K) * sizeof(*F));

	for (uint32_t i = 0; i < N_elems; i++)
		assemble_elem(i, K, F, part, graph, material,
			      enable_self_weight, gravity,
			      analysis2D, params2D);
	
	return 0;
}

static void assemble_elem(uint32_t elem_id,
			  vcn_sparse_t *K, double *F,
			  const nb_partition_t *const part,
			  const nb_graph_t *graph,
			  const nb_material_t *material,
			  bool enable_self_weight,
			  double gravity[2],
			  nb_analysis2D_t analysis2D,
			  nb_analysis2D_params *params2D)
{
	uint16_t N = 1 + graph->N_adj[elem_id];
	uint32_t memsize = (4 * N) * sizeof(double) +
		N * sizeof(uint32_t);
	if (ENABLE_LEAST_SQUARES)
		memsize += POW2(2 * N) * sizeof(double);

	char *memblock = NB_SOFT_MALLOC(memsize);
	uint32_t *ngb = (void*) memblock;
	double *Ke = (void*) (memblock + N * sizeof(uint32_t));

	integrate_elem_force(part, material, enable_self_weight,
			     gravity, elem_id, F);

	integrate_elem(elem_id, Ke, F, part, graph, material,
		       analysis2D, params2D);
	
	if (ENABLE_LEAST_SQUARES) {
		double *KTK = (void*) (memblock + N * sizeof(uint32_t) +
				       4 * N * sizeof(double));
		get_KTKe(N, Ke, KTK);
		add_KTKe_to_K(graph, elem_id, KTK, K);
	} else {
		add_Ke_to_K(graph, elem_id, Ke, K);
	}

	NB_SOFT_FREE(memsize, memblock);
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

static void integrate_elem(uint32_t elem_id,
			   double *Ke, double *F,
			   const nb_partition_t *const part,
			   const nb_graph_t *graph,
			   const nb_material_t *material,
			   nb_analysis2D_t analysis2D,
			   nb_analysis2D_params *params2D)
{
	uint32_t N = 1 + graph->N_adj[elem_id];
	memset(Ke, 0, 4 * N * sizeof(double));
	uint16_t N_faces = nb_partition_elem_get_N_adj(part, elem_id);
	for (uint16_t j = 0; j < N_faces; j++) {
		integrate_face(elem_id, j, part, graph, material,
			       analysis2D, params2D, Ke);
	}
}

static void integrate_face(uint32_t elem_id, 
			   uint16_t face_id,
			   const nb_partition_t *const part,
			   const nb_graph_t *graph,
			   const nb_material_t *material,
			   nb_analysis2D_t analysis2D,
			   nb_analysis2D_params *params2D,
			   double *Ke)
{
	if (nb_partition_elem_has_ngb(part, elem_id, face_id)) {
		double D[4];
		nb_pde_get_constitutive_matrix(D, material, analysis2D);

		uint16_t N_ngb = face_get_N_ngb(part, elem_id, face_id);

		uint32_t memsize = 4 * N_ngb * sizeof(double) +
			N_ngb * sizeof(uint32_t);
		char* memblock = NB_SOFT_MALLOC(memsize);
		uint32_t *ngb = (void*) memblock;
		double *Kf = (void*) (memblock + N_ngb * sizeof(uint32_t));

		face_get_neighbourhood(part, elem_id, face_id, ngb);

		get_Kf(part, elem_id, face_id, D, N_ngb, ngb,
		       QUADRATURE_POINTS, params2D, Kf);

		add_Kf_to_Ke(Ke, Kf, graph, elem_id, N_ngb, ngb);

		NB_SOFT_FREE(memsize, memblock);
	} else {
		/* BOUNDARY CONDITION */
	}
}

static uint16_t face_get_N_ngb(const nb_partition_t *part,
			       uint32_t elem_id, uint16_t face_id)
{	
	uint16_t N = 2;
	if (ALL_NEIGHBOURS) {
		N += get_N_ngb_around_right_vtx(part, elem_id, face_id);
		uint32_t ngb_id = 
			nb_partition_elem_get_ngb(part, elem_id, face_id);
		uint16_t aux = nb_partition_elem_ngb_get_face(part, ngb_id,
							      elem_id);
		N += get_N_ngb_around_right_vtx(part, ngb_id, aux);
	}
	return N;
}

static uint16_t get_N_ngb_around_right_vtx(const nb_partition_t *part,
					   uint32_t elem_id,
					   uint16_t face_id)
{
	uint32_t N_elems = nb_partition_get_N_elems(part);

	uint32_t front_ngb_id = nb_partition_elem_get_ngb(part, elem_id,
							  face_id);
	uint32_t nid_prev = elem_id;
	uint32_t nid = nb_partition_elem_face_get_right_ngb(part, elem_id,
							    face_id);
	uint16_t N = 0;
	while (nid != front_ngb_id && nid < N_elems) {
		N += 1;
		uint16_t aux = nb_partition_elem_ngb_get_face(part, nid,
							      nid_prev);
		nid_prev = nid;
		nid = nb_partition_elem_face_get_right_ngb(part, nid, aux);
	}
	if (nid >= N_elems && front_ngb_id < N_elems) {
		nid_prev = front_ngb_id;
		uint16_t aux = nb_partition_elem_ngb_get_face(part, front_ngb_id,
							      elem_id);
		nid = nb_partition_elem_face_get_left_ngb(part, front_ngb_id,
							  aux);
		while (nid < N_elems) {
			N += 1;
			aux = nb_partition_elem_ngb_get_face(part, nid,
							     nid_prev);
			nid_prev = nid;
			nid = nb_partition_elem_face_get_left_ngb(part, nid,
								  aux);
		}
	}
	return N;
}

static void face_get_neighbourhood(const nb_partition_t *part,
				   uint32_t elem_id, uint16_t face_id,
				   uint32_t *ngb)
{
	ngb[0] = elem_id;
	ngb[1] = nb_partition_elem_get_ngb(part, elem_id, face_id);
	uint16_t id = 2;

	if (ALL_NEIGHBOURS) {
		id = get_ngb_around_right_vtx(part, ngb, id, ngb[0], face_id);
		uint16_t aux = nb_partition_elem_ngb_get_face(part, ngb[1],
							      ngb[0]);
		get_ngb_around_right_vtx(part, ngb, id, ngb[1], aux);
	}
}

static uint16_t get_ngb_around_right_vtx(const nb_partition_t *part,
					 uint32_t *ngb, uint16_t current_id,
					 uint32_t elem_id,
					 uint16_t face_id)
{
	uint32_t N_elems = nb_partition_get_N_elems(part);

	uint32_t front_ngb_id = nb_partition_elem_get_ngb(part, elem_id,
							  face_id);
	uint32_t nid_prev = elem_id;
	uint32_t nid = nb_partition_elem_face_get_right_ngb(part, elem_id,
							    face_id);
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
		nid = nb_partition_elem_face_get_left_ngb(part, front_ngb_id,
							  aux);
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
		   nb_analysis2D_params *params2D,
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
		double factor = wqp[q] * params2D->thickness;
		for (uint16_t i = 0; i < N_ngb; i++) {
			double Kfi[4];
			get_Kf_nodal_contribution(part, D, nf, i,
						  grad_phi, Kfi);
			Kf[i * 2] += factor * Kfi[0];
			Kf[i*2+1] += factor * Kfi[1];
			Kf[2 * N_ngb + i * 2] += factor * Kfi[2];
			Kf[2 * N_ngb + i*2+1] += factor * Kfi[3];
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

	nb_nonpolynomial_eval_grad(N_ngb, 2, ni, NULL, x, 0, grad_phi);	

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

static void add_Kf_to_Ke(double *Ke, const double *Kf,
			 const nb_graph_t *graph,  uint32_t elem_id,
			 uint16_t N_ngb, const uint32_t *ngb)
{
	uint16_t N = 1 + graph->N_adj[elem_id];
	for (uint32_t m = 0; m < N_ngb; m++) {
		uint32_t imsh = ngb[m];
		uint16_t i = get_id_local_from_global(graph, elem_id, imsh);
		Ke[i * 2] += Kf[2 * m];
		Ke[i*2+1] += Kf[2*m+1];
		Ke[2 * N + (i * 2)] += Kf[2 * N_ngb + (2 * m)];
		Ke[2 * N + (i*2+1)] += Kf[2 * N_ngb + (2*m+1)];
	}
}

static uint16_t get_id_local_from_global(const nb_graph_t *graph,
					 uint32_t elem_id,
					 uint32_t global_id)
{
	uint16_t id = 0;
	if (elem_id != global_id) {
		uint16_t N_adj = graph->N_adj[elem_id];
		for (uint16_t i = 0; i < N_adj; i++) {
			if (global_id == graph->adj[elem_id][i]) {
				id = i + 1;
				break;
			}
		}
	}
	return id;
}

static double get_KTKe(uint16_t N, const double *Ke, double *KTK)
{
	for (uint16_t i = 0; i < 2 * N; i++) {
		for (uint16_t j = 0; j < 2 * N; j++) {
			double dot = 0;
			for (uint8_t k = 0; k < 2; k++) {
				dot += Ke[k * 2 * N + i] *
					Ke[k * 2 * N + j];
			}
			KTK[i * 2 * N + j] = dot;
		}
	}
}

static void add_KTKe_to_K(const nb_graph_t *graph, uint32_t elem_id,
			  const double *KTK, vcn_sparse_t *K)
{
	uint16_t N = 1 + graph->N_adj[elem_id];
	uint16_t size = 2 * N;
	for (uint32_t m = 0; m < N; m++) {
		uint32_t i = get_id_global_from_local(graph, elem_id, m);
		for (uint32_t n = 0; n < N; n++) {
			uint32_t j = get_id_global_from_local(graph, 
							      elem_id, n);
			vcn_sparse_add(K, i * 2, j * 2,
				       KTK[(2 * m)*size + (2 * n)]);
			vcn_sparse_add(K, i * 2, j*2+1,
				       KTK[(2 * m)*size + (2*n+1)]);
			vcn_sparse_add(K, i*2+1, j * 2,
				       KTK[(2*m+1)*size + (2 * n)]);
			vcn_sparse_add(K, i*2+1, j*2+1,
				       KTK[(2*m+1)*size + (2*n+1)]);
		}
	}
}

static uint32_t get_id_global_from_local(const nb_graph_t *graph,
					 uint32_t elem_id,
					 uint16_t local_id)
{
	uint32_t id = elem_id;
	if (0 < local_id)
		id = graph->adj[elem_id][local_id - 1];
	return id;
}

static void add_Ke_to_K(const nb_graph_t *graph,  uint32_t elem_id,
			const double *Ke, vcn_sparse_t *K)
{
	uint16_t N = 1 + graph->N_adj[elem_id];
	uint32_t i = elem_id;
	uint16_t size = 2 * N;
	for (uint32_t m = 0; m < N; m++) {
		uint32_t j = get_id_global_from_local(graph, elem_id, m);
		vcn_sparse_add(K, i * 2, j * 2, Ke[m * 2]);
		vcn_sparse_add(K, i * 2, j*2+1, Ke[m*2+1]);
		vcn_sparse_add(K, i*2+1, j * 2, Ke[size + m * 2]);
		vcn_sparse_add(K, i*2+1, j*2+1, Ke[size + m*2+1]);
	}
}

static int solver(const vcn_sparse_t *const A,
		  const double *const b, double* x)
{
	int status;
	if (ENABLE_LEAST_SQUARES) {
		uint32_t N = vcn_sparse_get_size(A);
		status = vcn_sparse_solve_CG_precond_Jacobi(A, b, x, N,
							    1e-8, NULL,
							    NULL, 1);
	} else {
		status = vcn_sparse_solve_using_LU(A, b, x, 1);
	}
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
