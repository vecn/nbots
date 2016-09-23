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
#define ENABLE_PAIR_WISE false
#define ENABLE_LEAST_SQUARES true


#define POW2(a) ((a)*(a))

static int assemble_system(vcn_sparse_t *K, double *F,
			   const nb_partition_t *const part,
			   const nb_graph_t *face_elems_conn,
			   const nb_material_t *material,
			   bool enable_self_weight,
			   double gravity[2],
			   nb_analysis2D_t analysis2D,
			   nb_analysis2D_params *params2D,
			   uint8_t N_qp);

static void integrate_elem_force(const nb_partition_t *part,
				 const nb_material_t *material,
				 bool enable_self_weight,
				 double gravity[2],
				 uint32_t elem_id,
				 double *F);
static void load_face_elems_conn(nb_graph_t *face_elems_conn,
				 const nb_partition_t *part);
static void load_elem_faces_conn(nb_graph_t *elem_faces_conn,
				 const nb_graph_t *face_elems_conn,
				 const nb_partition_t *part);
static uint32_t get_N_total_face_adj(const nb_partition_t *part);
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
					 uint32_t elem_id, uint16_t face_id);
static void assemble_face(uint16_t face_id,
			  vcn_sparse_t *K, double *F,
			  const nb_partition_t *const part,
			  const nb_graph_t *face_elems_conn,
			  const nb_material_t *material,
			  nb_analysis2D_t analysis2D,
			  nb_analysis2D_params *params2D,
			  uint8_t N_qp);
static void integrate_Kf(const nb_partition_t *const part,
			 const double D[4], uint16_t N,
			 const uint32_t *adj, uint8_t N_qp,
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
static void add_KTKf_to_K(uint16_t N, uint32_t *adj,
			  const double *Kf, vcn_sparse_t *K);
static double get_KTKf(uint16_t N, const double *Kf, double *KTK);
static void add_Kf_to_K(uint16_t N, uint32_t *adj,
			const double *Kf, vcn_sparse_t *K);
static int solver(const vcn_sparse_t *const A,
		  const double *const b, double* x);

static void vector_permutation(uint32_t N, const double *v,
			       const uint32_t *perm, double *vp);
static void compute_strain(double *strain,
			   const nb_graph_t *face_elems_conn,
			   const nb_graph_t *elem_faces_conn,
			   const nb_partition_t *const part,
			   double *disp, uint8_t N_qp);
static void get_face_strain(uint32_t face_id,
			    const nb_graph_t *face_elems_conn,
			    const nb_partition_t *const part,
			    double *disp, uint8_t N_qp,
			    double *fstrain);
static void get_elem_strain(uint32_t elem_id, const double *fstrain,
			    const nb_graph_t *elem_faces_conn,
			    const nb_partition_t *const part,
			    double *disp, uint8_t N_qp, double *strain);

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
	uint8_t N_qp = QUADRATURE_POINTS;
	uint32_t N_elems = nb_partition_get_N_elems(part);
	uint32_t graph_size = nb_graph_get_memsize();
	uint32_t system_size = 2 * N_elems * sizeof(double);
	uint32_t memsize = 3 * graph_size + system_size;
	char *memblock = NB_SOFT_MALLOC(memsize);

	double* F = (void*) memblock;
	nb_graph_t *graph = (void*) (memblock + system_size);
	nb_graph_t *face_elems_conn = (void*)
		(memblock + graph_size + system_size);
	nb_graph_t *elem_faces_conn = (void*)
		(memblock + 2 * graph_size + system_size);

	nb_graph_init(graph);
	nb_partition_load_graph(part, graph, NB_ELEMS_LINKED_BY_NODES);
	nb_graph_extend_adj(graph, 1);
	vcn_sparse_t *K = vcn_sparse_create(graph, NULL, 2);
	nb_graph_finish(graph);

	memset(F, 0, system_size);

	nb_graph_init(face_elems_conn);
	load_face_elems_conn(face_elems_conn, part);

	nb_graph_init(elem_faces_conn);
	load_elem_faces_conn(elem_faces_conn, face_elems_conn, part);

	int status_assemble = assemble_system(K, F, part, face_elems_conn,
					      material, enable_self_weight,
					      gravity, analysis2D, params2D,
					      N_qp);


	if (0 != status_assemble) {
		status = 1;
		goto CLEANUP_LINEAR_SYSTEM;
	}

	nb_cvfa_set_bconditions(part, material, analysis2D,
				K, F, bcond, 1.0);

	int solver_status = solver(K, F, displacement);
	if (0 != solver_status) {
		status = 2;
		goto CLEANUP_LINEAR_SYSTEM;
	}
	
	compute_strain(strain, face_elems_conn, elem_faces_conn,
		       part, displacement, N_qp);

	status = 0;
CLEANUP_LINEAR_SYSTEM:
	vcn_sparse_destroy(K);
	nb_graph_finish(face_elems_conn);
	nb_graph_finish(elem_faces_conn);
	NB_SOFT_FREE(memsize, memblock);
	return status;
}

static int assemble_system(vcn_sparse_t *K, double *F,
			   const nb_partition_t *const part,
			   const nb_graph_t *face_elems_conn,
			   const nb_material_t *material,
			   bool enable_self_weight,
			   double gravity[2],
			   nb_analysis2D_t analysis2D,
			   nb_analysis2D_params *params2D,
			   uint8_t N_qp)
{
	vcn_sparse_reset(K);
	memset(F, 0, vcn_sparse_get_size(K) * sizeof(*F));

	uint32_t N_elems = nb_partition_get_N_elems(part);
	for (uint32_t i = 0; i < N_elems; i++)
		integrate_elem_force(part, material, enable_self_weight,
				     gravity, i, F);

	for (uint32_t i = 0; i < face_elems_conn->N; i++)
		assemble_face(i, K, F, part, face_elems_conn, material,
			      analysis2D, params2D, N_qp);

	return 0;
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

static void load_face_elems_conn(nb_graph_t *face_elems_conn,
				 const nb_partition_t *part)
{
	uint32_t N = nb_partition_get_N_edges(part);
	uint32_t N_total_adj = get_N_total_face_adj(part);
	face_elems_conn->N = N;
	uint32_t memsize = N * (sizeof(*(face_elems_conn->N_adj)) +
				sizeof(*(face_elems_conn->adj))) +
		N_total_adj * sizeof(**(face_elems_conn->adj));
	char *memblock = malloc(memsize);
	face_elems_conn->N_adj = (void*) memblock;
	face_elems_conn->adj = (void*)
		(memblock + N * sizeof(*(face_elems_conn->N_adj)));
	
	uint32_t id = 0;
	uint32_t N_elems = nb_partition_get_N_elems(part);
	memblock += N * (sizeof(*(face_elems_conn->N_adj)) +
			 sizeof(*(face_elems_conn->adj)));

	for (uint32_t i = 0; i < N_elems; i++) {
		uint16_t N_adj = nb_partition_elem_get_N_adj(part, i);
		for (uint32_t j = 0; j < N_adj; j++) {
			uint32_t ngb_id = nb_partition_elem_get_ngb(part, i, j);
			if (N_elems <= ngb_id) {
				uint16_t N_ngb = 1;
				face_elems_conn->N_adj[id] = N_ngb;
				face_elems_conn->adj[id] = (void*) memblock;
				memblock += N_ngb * sizeof(**(face_elems_conn->adj));
				face_elems_conn->adj[id][0] = i;
				id += 1;
			} else if (i < ngb_id) {
				uint16_t N_ngb = face_get_N_ngb(part, i, j);
				face_elems_conn->N_adj[id] = N_ngb;
				face_elems_conn->adj[id] = (void*) memblock;
				memblock += N_ngb * sizeof(**(face_elems_conn->adj));
				face_get_neighbourhood(part, i, j,
						       face_elems_conn->adj[id]);
				id += 1;
			}
		}
	}
}

static void load_elem_faces_conn(nb_graph_t *elem_faces_conn,
				 const nb_graph_t *face_elems_conn,
				 const nb_partition_t *part)
{
	uint32_t N = nb_partition_get_N_elems(part);
	uint32_t N_total_adj = 0;
	for (uint32_t i = 0; i < N; i++)
		N_total_adj += nb_partition_elem_get_N_adj(part, i);
	
	elem_faces_conn->N = N;
	uint32_t memsize = N * (sizeof(*(elem_faces_conn->N_adj)) +
				sizeof(*(elem_faces_conn->adj))) +
		N_total_adj * sizeof(**(elem_faces_conn->adj));
	char *memblock = malloc(memsize);
	elem_faces_conn->N_adj = (void*) memblock;
	elem_faces_conn->adj = (void*)
		(memblock + N * sizeof(*(elem_faces_conn->N_adj)));
	
	memblock += N * (sizeof(*(elem_faces_conn->N_adj)) +
			 sizeof(*(elem_faces_conn->adj)));
	for (uint32_t i = 0; i < N; i++) {
		uint16_t N_adj = nb_partition_elem_get_N_adj(part, i);
		elem_faces_conn->N_adj[i] = N_adj;
		elem_faces_conn->adj[i] = (void*) memblock;
		memblock += N_adj * sizeof(**(elem_faces_conn->adj));
	}

	memset(elem_faces_conn->N_adj, 0,
	       N * sizeof(*(elem_faces_conn->N_adj)));
	for (uint32_t i = 0; i < face_elems_conn->N; i++) {
		uint32_t elem_1 = face_elems_conn->adj[i][0];
		uint32_t elem_2 = face_elems_conn->adj[i][1];

		uint16_t face_1 =
			nb_partition_elem_ngb_get_face(part, elem_1, elem_2);
		uint16_t face_2 =
			nb_partition_elem_ngb_get_face(part, elem_2, elem_1);

		elem_faces_conn->adj[elem_1][face_1] = i;
		elem_faces_conn->adj[elem_2][face_2] = i;
	}	
}

static uint32_t get_N_total_face_adj(const nb_partition_t *part)
{
	uint32_t N = 0;
	uint32_t N_elems = nb_partition_get_N_elems(part);
	for (uint32_t i = 0; i < N_elems; i++) {
		uint16_t N_adj = nb_partition_elem_get_N_adj(part, i);
		for (uint32_t j = 0; j < N_adj; j++) {
			uint32_t ngb_id = nb_partition_elem_get_ngb(part,
								    i, j);
			if (N_elems <= ngb_id)
				N += 1;
			else if (i < ngb_id)
				N += face_get_N_ngb(part, i, j);
		}
	}
	return N;
}

static uint16_t face_get_N_ngb(const nb_partition_t *part,
			       uint32_t elem_id, uint16_t face_id)
{	
	uint16_t N = 2;
	if (!ENABLE_PAIR_WISE) {
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
	if (!ENABLE_PAIR_WISE) {
		uint16_t id = 2;
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

static void assemble_face(uint16_t face_id,
			  vcn_sparse_t *K, double *F,
			  const nb_partition_t *const part,
			  const nb_graph_t *face_elems_conn,
			  const nb_material_t *material,
			  nb_analysis2D_t analysis2D,
			  nb_analysis2D_params *params2D,
			  uint8_t N_qp)
{
	uint16_t N = face_elems_conn->N_adj[face_id];
	if (N > 1) {
		double D[4];
		nb_pde_get_constitutive_matrix(D, material, analysis2D);

		uint32_t *adj = face_elems_conn->adj[face_id];
		uint32_t memsize = 4 * N * sizeof(double);
		char* memblock = NB_SOFT_MALLOC(memsize);
		double *Kf = (void*) memblock;

		integrate_Kf(part, D, N, adj, N_qp, params2D, Kf);

		if (ENABLE_LEAST_SQUARES)
			add_KTKf_to_K(N, adj, Kf, K);
		else
			add_Kf_to_K(N, adj, Kf, K);

		NB_SOFT_FREE(memsize, memblock);
	}
}

static void integrate_Kf(const nb_partition_t *const part,
			 const double D[4], uint16_t N,
			 const uint32_t *adj, uint8_t N_qp,
			 nb_analysis2D_params *params2D,
			 double *Kf)
{
	uint32_t memsize = (3 * N_qp + 2 * N) * sizeof(double);
	char *memblock = NB_SOFT_MALLOC(memsize);
	double *wqp = (void*) memblock;
	double *xqp = (void*) (memblock + N_qp * sizeof(double));
	double *grad_phi = (void*) (memblock + 3 * N_qp * sizeof(double));

	double nf[2];
	uint16_t local_face_id = nb_partition_elem_ngb_get_face(part, adj[0],
								adj[1]);
	double lf = nb_partition_elem_face_get_normal(part, adj[0],
						      local_face_id, nf);

	get_quadrature_points(part, adj[0], local_face_id, lf,
			      N_qp, xqp, wqp);

	memset(Kf, 0, 4 * N * sizeof(double));
	for (uint8_t q = 0; q < N_qp; q++) {
		interpolators_eval_grad(part, N, adj, 
					&(xqp[q*2]), grad_phi);
		double factor = wqp[q] * params2D->thickness;
		for (uint16_t i = 0; i < N; i++) {
			double Kfi[4];
			get_Kf_nodal_contribution(part, D, nf, i,
						  grad_phi, Kfi);
			Kf[i * 2] += factor * Kfi[0];
			Kf[i*2+1] += factor * Kfi[1];
			Kf[2 * N + i * 2] += factor * Kfi[2];
			Kf[2 * N + i*2+1] += factor * Kfi[3];
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
		wqp[q] = lf * (glq.w[q] / 2.0);
		double w = (glq.x[q] + 1)/2.0;
		nb_partition_elem_face_get_midpoint(part, elem_id, face_id,
						    w, &(xqp[q*2]));
	}
}

static void interpolators_eval_grad(const nb_partition_t *part, uint8_t N_ngb,
				    const uint32_t *ngb, const double x[2],
				    double *grad_phi)
{	
	uint32_t memsize = 3 * N_ngb * sizeof(double);
	char *memblock = NB_SOFT_MALLOC(memsize);
	double *ni = (void*) memblock;
	double *ri = (void*) (memblock + 2 * N_ngb * sizeof(double));

	for (uint32_t i = 0; i < N_ngb; i++) {
		ni[i * 2] = nb_partition_elem_get_x(part, ngb[i]);
		ni[i*2+1] = nb_partition_elem_get_y(part, ngb[i]);
		ri[i] = nb_partition_elem_get_apotem(part, ngb[i]);
	}

	nb_nonpolynomial_eval_grad(N_ngb, 2, ni, ri, x, grad_phi);	

	NB_SOFT_FREE(memsize, memblock);
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

static void add_KTKf_to_K(uint16_t N, uint32_t *adj,
			  const double *Kf, vcn_sparse_t *K)
{
	uint32_t memsize = POW2(2 * N) * sizeof(double);
	char* memblock = NB_SOFT_MALLOC(memsize);
	double *KTK = (void*) memblock;
	get_KTKf(N, Kf, KTK);

	uint16_t size = 2 * N;
	for (uint32_t m = 0; m < N; m++) {
		uint32_t i = adj[m];
		for (uint32_t n = 0; n < N; n++) {
			uint32_t j = adj[n];
			vcn_sparse_add(K, i * 2, j * 2,
				       2 * KTK[(2 * m)*size + (2 * n)]);
			vcn_sparse_add(K, i * 2, j*2+1,
				       2 * KTK[(2 * m)*size + (2*n+1)]);
			vcn_sparse_add(K, i*2+1, j * 2,
				       2 * KTK[(2*m+1)*size + (2 * n)]);
			vcn_sparse_add(K, i*2+1, j*2+1,
				       2 * KTK[(2*m+1)*size + (2*n+1)]);
		}
	}
	NB_SOFT_FREE(memsize, memblock);
}

static double get_KTKf(uint16_t N, const double *Kf, double *KTK)
{
	for (uint16_t i = 0; i < 2 * N; i++) {
		for (uint16_t j = 0; j < 2 * N; j++) {
			double dot = 0;
			for (uint8_t k = 0; k < 2; k++) {
				dot += Kf[k * 2 * N + i] *
					Kf[k * 2 * N + j];
			}
			KTK[i * 2 * N + j] = dot;
		}
	}
}

static void add_Kf_to_K(uint16_t N, uint32_t *adj,
			const double *Kf, vcn_sparse_t *K)
{
	uint16_t size = 2 * N;
	
	uint32_t i = adj[0];
	for (uint32_t m = 0; m < N; m++) {
		uint32_t j = adj[m];
		vcn_sparse_add(K, i * 2, j * 2, -Kf[m * 2]);
		vcn_sparse_add(K, i * 2, j*2+1, -Kf[m*2+1]);
		vcn_sparse_add(K, i*2+1, j * 2, -Kf[size + m * 2]);
		vcn_sparse_add(K, i*2+1, j*2+1, -Kf[size + m*2+1]);
	}

	i = adj[1];
	for (uint32_t m = 0; m < N; m++) {
		uint32_t j = adj[m];
		vcn_sparse_add(K, i * 2, j * 2, Kf[m * 2]);
		vcn_sparse_add(K, i * 2, j*2+1, Kf[m*2+1]);
		vcn_sparse_add(K, i*2+1, j * 2, Kf[size + m * 2]);
		vcn_sparse_add(K, i*2+1, j*2+1, Kf[size + m*2+1]);
	}
}

static int solver(const vcn_sparse_t *const A,
		  const double *const b, double* x)
{
	uint32_t N = vcn_sparse_get_size(A);
	uint16_t grp_size = nb_graph_get_memsize();
	uint32_t memsize = grp_size +
		2 * N * (sizeof(uint32_t) + sizeof(double));
	char *memblock = NB_SOFT_MALLOC(memsize);
	nb_graph_t *graph = (void*) memblock;
	uint32_t *perm = (void*) (memblock + grp_size);
	uint32_t *iperm = (void*) (memblock + grp_size +
				   N * sizeof(uint32_t));
	double *br = (void*) (memblock + grp_size +
			      2 * N * sizeof(uint32_t));
	double *xr = (void*) (memblock + grp_size +
			      2 * N * sizeof(uint32_t) + N * sizeof(double));

	nb_graph_init(graph);
	nb_sparse_get_graph(A, graph);
	nb_graph_labeling(graph, perm, iperm, NB_LABELING_ND);
	nb_graph_finish(graph);

	vcn_sparse_t *Ar = vcn_sparse_create_permutation(A, perm, iperm);
	vector_permutation(N, b, perm, br);

	int status = vcn_sparse_solve_Cholesky(Ar, br, xr, 1);
	if (0 != status)
		status = vcn_sparse_solve_using_LU(Ar, br, xr, 1);

	vector_permutation(N, xr, iperm, x);
	
	NB_SOFT_FREE(memsize, memblock);
	return status;
}

static void vector_permutation(uint32_t N, const double *v,
			       const uint32_t *perm, double *vp)
{
	for (uint32_t i = 0; i < N; i++)
		vp[i] = v[perm[i]];
}

#include "nb/pde_bot/frechet_derivative.h"/* TEMPORAL */
static void compute_strain(double *strain,
			   const nb_graph_t *face_elems_conn,
			   const nb_graph_t *elem_faces_conn,
			   const nb_partition_t *const part,
			   double *disp, uint8_t N_qp)
{
	uint32_t N_elems = nb_partition_get_N_elems(part);
	
	nb_graph_t *graph = alloca(nb_graph_get_memsize());
	nb_graph_init(graph);
	nb_partition_load_graph(part, graph, NB_ELEMS_LINKED_BY_NODES);
	nb_graph_extend_adj(graph, 6);

	double x[2], u[2], xi[6000], ui[6000];
	for (uint32_t i = 0; i < graph->N; i++) {
		x[0] = nb_partition_elem_get_x(part, i);
		x[1] = nb_partition_elem_get_y(part, i);
		u[0] = disp[i * 2];
		u[1] = disp[i*2+1];
		uint16_t N = 0;
		uint16_t N_adj = graph->N_adj[i];
		for (uint16_t j = 0; j < N_adj; j++) {
			uint32_t ngb = graph->adj[i][j];
			xi[N * 2] = nb_partition_elem_get_x(part, ngb);
			xi[N*2+1] = nb_partition_elem_get_y(part, ngb);
			ui[N * 2] = disp[ngb * 2];
			ui[N*2+1] = disp[ngb*2+1];
			N += 1;
		}
		double Du[4];
		nb_pde_get_frechet_derivative(N, 2, 2, x, u, xi, ui, Du);
		strain[i * 3] = Du[0];
		strain[i*3+1] = Du[3];
		strain[i*3+2] = Du[1] + Du[2];
	}
	nb_graph_finish(graph);
	return;/* AQUI VOY: Use all sourrounding 
		  elements (interpolation partners) */

	uint32_t memsize = 3 * face_elems_conn->N * N_qp * sizeof(double);
	double *fstrain = NB_SOFT_MALLOC(memsize);

	for (uint32_t i = 0; i < face_elems_conn->N; i++)
		get_face_strain(i, face_elems_conn, part,
				disp, N_qp, fstrain);

	for (uint32_t i = 0; i < N_elems; i++)
		get_elem_strain(i, fstrain, elem_faces_conn,
				part, disp, N_qp, strain);

	NB_SOFT_FREE(memsize, fstrain);
}

static void get_face_strain(uint32_t face_id,
			    const nb_graph_t *face_elems_conn,
			    const nb_partition_t *const part,
			    double *disp, uint8_t N_qp,
			    double *fstrain)
{
	uint16_t N = face_elems_conn->N_adj[face_id];
	if (N > 1) {
		uint32_t *adj = face_elems_conn->adj[face_id];
		uint32_t memsize = (3 * N_qp + 2 * N) * sizeof(double);
		char *memblock = NB_SOFT_MALLOC(memsize);
		double *wqp = (void*) memblock;
		double *xqp = (void*) (memblock + N_qp * sizeof(double));
		double *grad_phi = (void*) (memblock + 3 * N_qp * sizeof(double));

		uint16_t local_face_id = 
			nb_partition_elem_ngb_get_face(part, adj[0], adj[1]);
		double lf = nb_partition_elem_face_get_length(part, adj[0],
							      local_face_id);

		get_quadrature_points(part, adj[0], local_face_id,
				      lf, N_qp, xqp, wqp);

		for (uint8_t q = 0; q < N_qp; q++) {
			interpolators_eval_grad(part, N, adj,
						&(xqp[q*2]), grad_phi);
			uint32_t id = face_id * N_qp + q;
			memset(&(fstrain[id*3]), 0, 3 * sizeof(*fstrain));
			for (uint16_t i = 0; i < N; i++) {
				uint32_t elem_id = adj[i];
				double u = disp[elem_id * 2];
				double v = disp[elem_id*2+1];
				double dphi_dx = grad_phi[i * 2];
				double dphi_dy = grad_phi[i*2+1];
				fstrain[id * 3] += dphi_dx * u;
				fstrain[id*3+1] += dphi_dy * v;
				fstrain[id*3+2] += dphi_dy * u + dphi_dx * v;
			}
		}
		NB_SOFT_FREE(memsize, memblock);
	} else {
		for (uint8_t q = 0; q < N_qp; q++) {
			uint32_t id = face_id * N_qp + q;
			memset(&(fstrain[id*3]), 0, 3 * sizeof(*fstrain));
		}
		;/* Boundary condition TEMPORAL */
	}
}

static void get_elem_strain(uint32_t elem_id, const double *fstrain,
			    const nb_graph_t *elem_faces_conn,
			    const nb_partition_t *const part,
			    double *disp, uint8_t N_qp, double *strain)
{
	uint16_t N_adj = nb_partition_elem_get_N_adj(part, elem_id);
	uint16_t N = N_qp * N_adj;
	
	uint32_t memsize = (N_qp + 3 * N) * sizeof(double);
	char *memblock = NB_SOFT_MALLOC(memsize);
	double *phi = (void*) memblock;
	double *xi = (void*) (memblock + N * sizeof(double));
	double *wqp = (void*) (memblock + 3 * N * sizeof(double));

	for (uint16_t i = 0; i < N_adj; i++) {
		double lf = nb_partition_elem_face_get_length(part,
							      elem_id, i);
		get_quadrature_points(part, elem_id, i, lf, N_qp,
				      &(xi[i * N_qp * 2]), wqp);
	}

	double x[2];
	x[0] = nb_partition_elem_get_x(part, elem_id);
	x[1] = nb_partition_elem_get_y(part, elem_id);
	nb_nonpolynomial_eval(N, 2, xi, NULL, x, phi);

	memset(&(strain[elem_id * 3]), 0, 3 * sizeof(*strain));
	for (uint16_t i = 0; i < N_adj; i++) {
		uint32_t face_id = elem_faces_conn->adj[elem_id][i];
		for (uint8_t q = 0; q < N_qp; q++) {
			uint16_t id = i * N_qp + q;
			uint32_t fs_id = face_id * N_qp + q;
			double factor = phi[id];
			strain[elem_id * 3] += factor * fstrain[fs_id * 3];
			strain[elem_id*3+1] += factor * fstrain[fs_id*3+1];
			strain[elem_id*3+2] += factor * fstrain[fs_id*3+2];
		}
	}
	NB_SOFT_FREE(memsize, memblock);
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
