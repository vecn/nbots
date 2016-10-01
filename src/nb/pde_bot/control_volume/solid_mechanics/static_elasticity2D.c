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

/* DEFINE RBF FOR INTERPOLATOR: x * log(x+1) */
/* DEFINE QUADRATURE POINTS 3 */

#define NGB_LAYERS 0 /* Negative: Pair-wise */

#define POW2(a) ((a)*(a))

static uint32_t get_cvfa_memsize(uint32_t N_elems);
static void distribute_cvfa_memory(char *memblock,
				   uint32_t N_elems, double **F,
				   nb_graph_t **face_elems_conn,
				   nb_graph_t **ngb_graph);
static void load_neighbourhood_graph(const nb_partition_t *part,
				     nb_graph_t *ngb_graph);
static void init_global_matrix(const nb_graph_t *ngb_graph,
			       vcn_sparse_t **K);
static void load_face_elems_conn(const nb_partition_t *part,
				 const nb_graph_t *ngb_graph,
				 nb_graph_t *face_elems_conn);
static uint32_t get_N_total_face_adj(const nb_partition_t *part,
				     const nb_graph_t *ngb_graph);
static uint16_t face_get_N_ngb(const nb_partition_t *part,
			       const nb_graph_t *ngb_graph,
			       uint32_t elem_id1, uint32_t elem_id2);
static uint16_t get_N_ngb_around_right_vtx(const nb_partition_t *part,
					   uint32_t elem_id1,
					   uint32_t elem_id2);
static char *set_elemental_faces(nb_graph_t *face_elems_conn,
				 const nb_partition_t *part,
				 const nb_graph_t *ngb_graph,
				 char *memblock, uint32_t elem_id);
static char *set_boundary_face(nb_graph_t *face_elems_conn,
			       const nb_partition_t *part, 
			       char *memblock, uint32_t face_id,
			       uint32_t elem_id);
static char *set_internal_face(nb_graph_t *face_elems_conn,
			       const nb_partition_t *part,
			       const nb_graph_t *ngb_graph,
			       char *memblock, uint32_t face_id,
			       uint32_t elem_id1, uint32_t elem_id2);
static void set_pairwise_ngb(uint32_t *adj, uint32_t elem_id1,
			     uint32_t elem_id2);
static void set_local_ngb(uint32_t *adj, const nb_partition_t *part,
			  uint32_t elem_id1, uint32_t elem_id2);
static uint16_t get_ngb_around_right_vtx(const nb_partition_t *part,
					 uint32_t *ngb, uint16_t current_id,
					 uint32_t elem_id1, uint32_t elem_id2);
static void set_extended_ngb(uint32_t *adj, const nb_graph_t *ngb_graph,
			     uint32_t elem_id1, uint32_t elem_id2);
static void assemble_global_forces(double *F,
				   const nb_partition_t *const part,
				   const nb_material_t *material,
				   bool enable_self_weight,
				   double gravity[2]);
static void integrate_elem_force(const nb_partition_t *part,
				 const nb_material_t *material,
				 bool enable_self_weight,
				 double gravity[2],
				 uint32_t elem_id,
				 double *F);
static void assemble_global_stiffness(vcn_sparse_t *K,
				      const nb_partition_t *const part,
				      const nb_graph_t *face_elems_conn,
				      const nb_material_t *material,
				      nb_analysis2D_t analysis2D,
				      nb_analysis2D_params *params2D,
				      uint8_t N_qp);
static void assemble_face(uint32_t face_id, vcn_sparse_t *K,
			  const nb_partition_t *const part,
			  const nb_graph_t *face_elems_conn,
			  const nb_material_t *material,
			  nb_analysis2D_t analysis2D,
			  nb_analysis2D_params *params2D,
			  uint8_t N_qp);
static void assemble_internal_face(uint32_t face_id, vcn_sparse_t *K,
				   const nb_partition_t *const part,
				   const nb_graph_t *face_elems_conn,
				   const double D[4],
				   nb_analysis2D_params *params2D,
				   uint8_t N_qp);
static void integrate_Kf(const nb_partition_t *const part, const double D[4],
			 uint32_t face_id, uint16_t N,  const uint32_t *adj,
			 uint8_t N_qp, nb_analysis2D_params *params2D,
			 double *Kf);
static void get_quadrature_points(const nb_partition_t *part, uint32_t face_id,
				  double lf, uint8_t N_qp, double *xqp,
				  double *wqp);
static void interpolators_eval_grad(const nb_partition_t *part,
				    uint32_t face_id, double nf[2],
				    uint8_t N_ngb, const uint32_t *ngb,
				    const double x[2], double *grad_phi);
static void get_Kf_nodal_contribution(const nb_partition_t *part,
				      const double D[4], const double nf[2],
				      uint16_t i, const double *grad_phi,
				      double Kfi[4]);
static void add_Kf_to_K(uint16_t N, const uint32_t *adj,
			const double *Kf, vcn_sparse_t *K);
static int solver(const vcn_sparse_t *const A,
		  const double *const b, double* x);
static void get_permutation(const vcn_sparse_t *const A,
			    uint32_t *perm, uint32_t *iperm);
static void vector_permutation(uint32_t N, const double *v,
			       const uint32_t *perm, double *vp);
static void compute_strain(double *strain, char *boundary_mask,
			   const nb_graph_t *face_elems_conn,
			   const nb_partition_t *const part,
			   const nb_bcond_t *const bcond,
			   const double *disp, uint8_t N_qp);
static void get_face_strain(uint32_t face_id,
			    const nb_graph_t *face_elems_conn,
			    const nb_partition_t *const part,
			    const nb_bcond_t *const bcond,
			    const double *disp, uint8_t N_qp,
			    double *strain,
			    char *boundary_mask);
static void get_internal_face_strain(uint32_t face_id,
				     const nb_graph_t *face_elems_conn,
				     const nb_partition_t *const part,
				     const double *disp, uint8_t N_qp,
				     double length, double nf[2],
				     double *strain);
static void get_boundary_face_strain(uint32_t face_id,
				     const nb_graph_t *face_elems_conn,
				     const nb_partition_t *const part,
				     const nb_bcond_t *const bcond,
				     const double *disp, uint8_t N_qp,
				     double length, double *strain);

int nb_cvfa_compute_2D_Solid_Mechanics
			(const nb_partition_t *const part,
			 const nb_material_t *const material,
			 const nb_bcond_t *const bcond,
			 bool enable_self_weight, double gravity[2],
			 nb_analysis2D_t analysis2D,
			 nb_analysis2D_params *params2D,
			 uint8_t N_quadrature_points,
			 double *displacement, /* Output */
			 double *strain,       /* Output */
			 char *boundary_mask   /* Output */)
{
	int status;
	uint32_t N_elems = nb_partition_get_N_elems(part);
	uint32_t memsize = get_cvfa_memsize(N_elems);
	char *memblock = NB_SOFT_MALLOC(memsize);
	double *F;
	nb_graph_t *face_elems_conn;
	nb_graph_t *ngb_graph;
	distribute_cvfa_memory(memblock, N_elems, &F,
			       &face_elems_conn, &ngb_graph);

	nb_graph_init(ngb_graph);
	load_neighbourhood_graph(part, ngb_graph);

	vcn_sparse_t *K;
	init_global_matrix(ngb_graph, &K);

	nb_graph_init(face_elems_conn);
	load_face_elems_conn(part, ngb_graph, face_elems_conn);
	nb_graph_finish(ngb_graph);

	assemble_global_forces(F, part, material, enable_self_weight,
			       gravity);

	assemble_global_stiffness(K, part, face_elems_conn, material,
				  analysis2D, params2D,
				  N_quadrature_points);
	
	nb_cvfa_set_bconditions(part, material, analysis2D, 
				K, F, bcond, 1.0);

	int solver_status = solver(K, F, displacement);
	if (0 != solver_status) {
		status = 1;
		goto CLEANUP_LINEAR_SYSTEM;
	}

	compute_strain(strain, boundary_mask, face_elems_conn, part,
		       bcond, displacement, N_quadrature_points);

	status = 0;
CLEANUP_LINEAR_SYSTEM:
	vcn_sparse_destroy(K);
	nb_graph_finish(face_elems_conn);
	NB_SOFT_FREE(memsize, memblock);
	return status;
}

static uint32_t get_cvfa_memsize(uint32_t N_elems)
{
	uint32_t graph_size = nb_graph_get_memsize();
	uint32_t system_size = 2 * N_elems * sizeof(double);
	return 2 * graph_size + system_size;
}

static void distribute_cvfa_memory(char *memblock,
				   uint32_t N_elems, double **F,
				   nb_graph_t **face_elems_conn,
				   nb_graph_t **ngb_graph)
{
	uint32_t graph_size = nb_graph_get_memsize();
	uint32_t system_size = 2 * N_elems * sizeof(double);
	*F = (void*) memblock;
	*face_elems_conn = (void*) (memblock + system_size);
	*ngb_graph = (void*) (memblock + system_size + graph_size);
}

static void load_neighbourhood_graph(const nb_partition_t *part,
				     nb_graph_t *ngb_graph)
{
	if (0 > NGB_LAYERS) {
		nb_partition_load_graph(part, ngb_graph, NB_ELEMS_LINKED_BY_EDGES);
	} else if (0 == NGB_LAYERS) {
		nb_partition_load_graph(part, ngb_graph, NB_ELEMS_LINKED_BY_NODES);
	} else {
		nb_partition_load_graph(part, ngb_graph, NB_ELEMS_LINKED_BY_NODES);
		nb_graph_extend_adj(ngb_graph, NGB_LAYERS - 1);
	}
}

static void init_global_matrix(const nb_graph_t *ngb_graph,
			       vcn_sparse_t **K)
{
	uint32_t memsize = nb_graph_get_memsize();
	nb_graph_t *graph = NB_SOFT_MALLOC(memsize);
	nb_graph_copy(graph, ngb_graph);
	nb_graph_extend_adj(graph, 1);

	*K = vcn_sparse_create(graph, NULL, 2);

	nb_graph_finish(graph);
	NB_SOFT_FREE(memsize, graph);
}

static void load_face_elems_conn(const nb_partition_t *part,
				 const nb_graph_t *ngb_graph,
				 nb_graph_t *face_elems_conn)
{
	uint32_t N = nb_partition_get_N_edges(part);
	uint32_t N_total_adj = get_N_total_face_adj(part, ngb_graph);
	face_elems_conn->N = N;
	uint32_t memsize = N * (sizeof(*(face_elems_conn->N_adj)) +
				sizeof(*(face_elems_conn->adj))) +
		N_total_adj * sizeof(**(face_elems_conn->adj));
	char *memblock = malloc(memsize);
	face_elems_conn->N_adj = (void*) memblock;
	face_elems_conn->adj = (void*)
		(memblock + N * sizeof(*(face_elems_conn->N_adj)));
	
	uint32_t N_elems = nb_partition_get_N_elems(part);
	memblock += N * (sizeof(*(face_elems_conn->N_adj)) +
			 sizeof(*(face_elems_conn->adj)));

	for (uint32_t i = 0; i < N_elems; i++)		
		memblock = set_elemental_faces(face_elems_conn, part,
					       ngb_graph, memblock, i);
}

static uint32_t get_N_total_face_adj(const nb_partition_t *part,
				     const nb_graph_t *ngb_graph)
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
				N += face_get_N_ngb(part, ngb_graph,
						    i, ngb_id);
		}
	}
	return N;
}

static uint16_t face_get_N_ngb(const nb_partition_t *part,
			       const nb_graph_t *ngb_graph,
			       uint32_t elem_id1, uint32_t elem_id2)
{
	uint16_t N;
	if (0 > NGB_LAYERS) {
		N = 2;
	} else if (0 == NGB_LAYERS) {
		N = 2 + get_N_ngb_around_right_vtx(part, elem_id1, elem_id2) +
			get_N_ngb_around_right_vtx(part, elem_id2, elem_id1);
	} else {
		uint16_t N_repeated =
			nb_graph_find_N_intersected_adj(ngb_graph, elem_id1,
							elem_id2);
		N = ngb_graph->N_adj[elem_id1] +
			ngb_graph->N_adj[elem_id2] - N_repeated;
	}
	return N;
}

static uint16_t get_N_ngb_around_right_vtx(const nb_partition_t *part,
					   uint32_t elem_id1,
					   uint32_t elem_id2)
{
	uint32_t N_elems = nb_partition_get_N_elems(part);

	uint32_t nid_prev = elem_id1;
	uint16_t aux = nb_partition_elem_ngb_get_face(part, elem_id1,
						      elem_id2);
	uint32_t nid = nb_partition_elem_face_get_right_ngb(part, elem_id1,
							    aux);
	uint16_t N = 0;
	while (nid != elem_id2 && nid < N_elems) {
		N += 1;
		aux = nb_partition_elem_ngb_get_face(part, nid, nid_prev);
		nid_prev = nid;
		nid = nb_partition_elem_face_get_right_ngb(part, nid, aux);
	}
	if (nid >= N_elems && elem_id2 < N_elems) {
		nid_prev = elem_id2;
		aux = nb_partition_elem_ngb_get_face(part, elem_id2, elem_id1);
		nid = nb_partition_elem_face_get_left_ngb(part, elem_id2,
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

static char *set_elemental_faces(nb_graph_t *face_elems_conn,
				 const nb_partition_t *part,
				 const nb_graph_t *ngb_graph,
				 char *memblock, uint32_t elem_id)
{
	uint16_t N_elems = nb_partition_get_N_elems(part);
	uint16_t N_adj = nb_partition_elem_get_N_adj(part, elem_id);
	for (uint32_t j = 0; j < N_adj; j++) {
		uint32_t ngb_id = nb_partition_elem_get_ngb(part, elem_id, j);
		if (N_elems <= ngb_id) {
			uint32_t face_id = 
				nb_partition_elem_find_edge(part, elem_id, j);
			set_boundary_face(face_elems_conn, part,
					  memblock, face_id, elem_id);
		} else if (elem_id < ngb_id) {
			uint32_t face_id =
				nb_partition_elem_find_edge(part, elem_id, j);
			memblock = set_internal_face(face_elems_conn, part,
						     ngb_graph, memblock,
						     face_id, elem_id, ngb_id);
		}
	}
	return memblock;
}

static char *set_boundary_face(nb_graph_t *face_elems_conn,
			       const nb_partition_t *part, 
			       char *memblock, uint32_t face_id,
			       uint32_t elem_id)
{
	face_elems_conn->N_adj[face_id] = 1;
	face_elems_conn->adj[face_id] = (void*) memblock;
	memblock += sizeof(**(face_elems_conn->adj));
	face_elems_conn->adj[face_id][0] = elem_id;
	return memblock;
}

static char *set_internal_face(nb_graph_t *face_elems_conn,
			       const nb_partition_t *part,
			       const nb_graph_t *ngb_graph,
			       char *memblock, uint32_t face_id,
			       uint32_t elem_id1, uint32_t elem_id2)
{
       	uint16_t N_ngb = face_get_N_ngb(part, ngb_graph, elem_id1, elem_id2);
	face_elems_conn->N_adj[face_id] = N_ngb;
	face_elems_conn->adj[face_id] = (void*) memblock;
	memblock += N_ngb * sizeof(**(face_elems_conn->adj));

	if (0 > NGB_LAYERS)
		set_pairwise_ngb(face_elems_conn->adj[face_id],
				 elem_id1, elem_id2);
	else if (0 == NGB_LAYERS)
		set_local_ngb(face_elems_conn->adj[face_id],
			      part, elem_id1, elem_id2);
	else
		set_extended_ngb(face_elems_conn->adj[face_id],
				 ngb_graph, elem_id1, elem_id2);

	return memblock;

}

static void set_pairwise_ngb(uint32_t *adj, uint32_t elem_id1,
			     uint32_t elem_id2)
{
	adj[0] = elem_id1;
	adj[1] = elem_id2;
}

static void set_local_ngb(uint32_t *adj, const nb_partition_t *part,
			  uint32_t elem_id1, uint32_t elem_id2)
{
	adj[0] = elem_id1;
	adj[1] = elem_id2;

	uint16_t id = get_ngb_around_right_vtx(part, adj, 2,
					       elem_id1, elem_id2);
	get_ngb_around_right_vtx(part, adj, id, elem_id2, elem_id1);
}

static uint16_t get_ngb_around_right_vtx(const nb_partition_t *part,
					 uint32_t *ngb, uint16_t current_id,
					 uint32_t elem_id1, uint32_t elem_id2)
{
	uint32_t N_elems = nb_partition_get_N_elems(part);

	uint32_t nid_prev = elem_id1;
	uint16_t aux = nb_partition_elem_ngb_get_face(part, elem_id1, elem_id2);
	uint32_t nid = nb_partition_elem_face_get_right_ngb(part, elem_id1,
							    aux);
	while (nid != elem_id2 && nid < N_elems) {
		ngb[current_id] = nid;
		current_id += 1;
		aux = nb_partition_elem_ngb_get_face(part, nid, nid_prev);
		nid_prev = nid;
		nid = nb_partition_elem_face_get_right_ngb(part, nid, aux);
	}
	if (nid >= N_elems && elem_id2 < N_elems) {
		nid_prev = elem_id2;
		aux = nb_partition_elem_ngb_get_face(part, elem_id2, elem_id1);
		nid = nb_partition_elem_face_get_left_ngb(part, elem_id2,
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

static void set_extended_ngb(uint32_t *adj, const nb_graph_t *ngb_graph,
			     uint32_t elem_id1, uint32_t elem_id2)
{
	adj[0] = elem_id1;
	adj[1] = elem_id2;
	uint16_t id = 2;
	for (uint16_t i = 0; i < ngb_graph->N_adj[elem_id1]; i++) {
		uint32_t ngb = ngb_graph->adj[elem_id1][i];
		if (ngb != elem_id2) {
			adj[id] = ngb;
			id += 1;
		}
	}
	for (uint16_t i = 0; i < ngb_graph->N_adj[elem_id2]; i++) {
		uint32_t ngb = ngb_graph->adj[elem_id2][i];
		bool included = false;
		for (uint16_t j = 0; j < id; j++) {	
			if (adj[j] == ngb) {
				included = true;
				break;
			}
		}
		if (!included) {
			adj[id] = ngb;
			id += 1;
		}
	}
}

static void assemble_global_forces(double *F,
				   const nb_partition_t *const part,
				   const nb_material_t *material,
				   bool enable_self_weight,
				   double gravity[2])
{
	uint32_t N_elems = nb_partition_get_N_elems(part);
	memset(F, 0, N_elems * 2 * sizeof(*F));
	for (uint32_t i = 0; i < N_elems; i++) {
		integrate_elem_force(part, material, enable_self_weight,
				     gravity, i, F);
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

static void assemble_global_stiffness(vcn_sparse_t *K,
				      const nb_partition_t *const part,
				      const nb_graph_t *face_elems_conn,
				      const nb_material_t *material,
				      nb_analysis2D_t analysis2D,
				      nb_analysis2D_params *params2D,
				      uint8_t N_qp)
{
	vcn_sparse_reset(K);
	uint32_t N_faces = nb_partition_get_N_edges(part);
	for (uint32_t i = 0; i < N_faces; i++) {
		assemble_face(i, K, part, face_elems_conn, material,
			      analysis2D, params2D, N_qp);
	}
}

static void assemble_face(uint32_t face_id, vcn_sparse_t *K,
			  const nb_partition_t *const part,
			  const nb_graph_t *face_elems_conn,
			  const nb_material_t *material,
			  nb_analysis2D_t analysis2D,
			  nb_analysis2D_params *params2D,
			  uint8_t N_qp)
{
	double D[4];
	nb_pde_get_constitutive_matrix(D, material, analysis2D);

	uint16_t N = face_elems_conn->N_adj[face_id];
	if (1 < N)
		assemble_internal_face(face_id, K, part, face_elems_conn,
				       D, params2D, N_qp);
}

static void assemble_internal_face(uint32_t face_id, vcn_sparse_t *K,
				   const nb_partition_t *const part,
				   const nb_graph_t *face_elems_conn,
				   const double D[4],
				   nb_analysis2D_params *params2D,
				   uint8_t N_qp)
{
	uint16_t N = face_elems_conn->N_adj[face_id];
	uint32_t *adj = face_elems_conn->adj[face_id];

	uint32_t memsize = 4 * N * sizeof(double);
	char* memblock = NB_SOFT_MALLOC(memsize);
	double *Kf = (void*) memblock;

	integrate_Kf(part, D, face_id, N, adj, N_qp, params2D, Kf);

	add_Kf_to_K(N, adj, Kf, K);

	NB_SOFT_FREE(memsize, memblock);
}

static void integrate_Kf(const nb_partition_t *const part, const double D[4],
			 uint32_t face_id, uint16_t N,  const uint32_t *adj,
			 uint8_t N_qp, nb_analysis2D_params *params2D,
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

	get_quadrature_points(part, face_id, lf, N_qp, xqp, wqp);

	memset(Kf, 0, 4 * N * sizeof(double));
	for (uint8_t q = 0; q < N_qp; q++) {
		interpolators_eval_grad(part, face_id, nf, N, adj, 
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

static void get_quadrature_points(const nb_partition_t *part, uint32_t face_id,
				  double lf, uint8_t N_qp, double *xqp,
				  double *wqp)
{
	nb_glquadrature_t glq;
	glq.x = alloca(N_qp * sizeof(*(glq.x)));
	glq.w = alloca(N_qp * sizeof(*(glq.w)));
	nb_glquadrature_load(&glq, N_qp);

	for (uint8_t q = 0; q < N_qp; q++) {
		wqp[q] = lf * (glq.w[q] / 2.0);
		double w = (glq.x[q] + 1)/2.0;
		nb_partition_edge_get_midpoint(part, face_id, w, &(xqp[q*2]));
	}
}

static void interpolators_eval_grad(const nb_partition_t *part,
				    uint32_t face_id, double nf[2],
				    uint8_t N_ngb, const uint32_t *ngb,
				    const double x[2], double *grad_phi)
{
	uint32_t memsize = 3 * N_ngb * sizeof(double);
	char *memblock = NB_SOFT_MALLOC(memsize);
	double *ni = (void*) memblock;
	double *ri = (void*) (memblock + 2 * N_ngb * sizeof(double));

	double r0 = nb_partition_elem_get_radius(part, ngb[0]);
	double r1 = nb_partition_elem_get_radius(part, ngb[1]);
	double r = vcn_math_harmonic_avg(r0, r1);
	for (uint16_t i = 0; i < N_ngb; i++) {
		ni[i * 2] = nb_partition_elem_get_x(part, ngb[i]);
		ni[i*2+1] = nb_partition_elem_get_y(part, ngb[i]);
		ri[i] = r;
	}
	ri[0] /= sqrt(2);
	ri[1] /= sqrt(2);

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

static void add_Kf_to_K(uint16_t N, const uint32_t *adj,
			const double *Kf, vcn_sparse_t *K)
{
	uint16_t size = 2 * N;
	uint32_t i = adj[0];
	uint32_t j = adj[1];
	for (uint32_t m = 0; m < N; m++) {
		uint32_t k = adj[m];
		vcn_sparse_add(K, i * 2, k * 2, -Kf[m * 2]);
		vcn_sparse_add(K, i * 2, k*2+1, -Kf[m*2+1]);
		vcn_sparse_add(K, i*2+1, k * 2, -Kf[size + m * 2]);
		vcn_sparse_add(K, i*2+1, k*2+1, -Kf[size + m*2+1]);

		vcn_sparse_add(K, j * 2, k * 2, Kf[m * 2]);
		vcn_sparse_add(K, j * 2, k*2+1, Kf[m*2+1]);
		vcn_sparse_add(K, j*2+1, k * 2, Kf[size + m * 2]);
		vcn_sparse_add(K, j*2+1, k*2+1, Kf[size + m*2+1]);
	}
}

static int solver(const vcn_sparse_t *const A,
		  const double *const b, double* x)
{
	uint32_t N = vcn_sparse_get_size(A);
	uint32_t memsize = 2 * N * (sizeof(uint32_t) + sizeof(double));
	char *memblock = NB_SOFT_MALLOC(memsize);
	uint32_t *perm = (void*) memblock;
	uint32_t *iperm = (void*) (memblock + N * sizeof(uint32_t));
	double *br = (void*) (memblock + 2 * N * sizeof(uint32_t));
	double *xr = (void*) (memblock + 2 * N * sizeof(uint32_t) +
			      N * sizeof(double));

	double usym = vcn_sparse_get_usym(A);/* TEMPORAL */
	printf("-- K usym: %e\n", usym);     /* TEMPORAL */

	get_permutation(A, perm, iperm);

	vcn_sparse_t *Ar = vcn_sparse_create_permutation(A, perm, iperm);
	vector_permutation(N, b, perm, br);

	int status = vcn_sparse_solve_using_LU(Ar, br, xr, 1);

	vector_permutation(N, xr, iperm, x);
	
	vcn_sparse_destroy(Ar);
	NB_SOFT_FREE(memsize, memblock);
	return status;
}

static void get_permutation(const vcn_sparse_t *const A,
			    uint32_t *perm, uint32_t *iperm)
{
	uint16_t memsize = nb_graph_get_memsize();
	nb_graph_t *graph = NB_SOFT_MALLOC(memsize);
	nb_graph_init(graph);
	nb_sparse_get_graph(A, graph);
	nb_graph_labeling(graph, perm, iperm, NB_LABELING_ND);
	nb_graph_finish(graph);

	NB_SOFT_FREE(memsize, graph);
}

static void vector_permutation(uint32_t N, const double *v,
			       const uint32_t *perm, double *vp)
{
	for (uint32_t i = 0; i < N; i++)
		vp[i] = v[perm[i]];
}

static void compute_strain(double *strain, char *boundary_mask,
			   const nb_graph_t *face_elems_conn,
			   const nb_partition_t *const part,
			   const nb_bcond_t *const bcond,
			   const double *disp, uint8_t N_qp)
{
	uint32_t N_faces = nb_partition_get_N_edges(part);

 	for (uint32_t i = 0; i < N_faces; i++) {
		get_face_strain(i, face_elems_conn, part, bcond,
				disp, N_qp, strain, boundary_mask);
	}
}

static void get_face_strain(uint32_t face_id,
			    const nb_graph_t *face_elems_conn,
			    const nb_partition_t *const part,
			    const nb_bcond_t *const bcond,
			    const double *disp, uint8_t N_qp,
			    double *strain,
			    char *boundary_mask)
{
	uint16_t N = face_elems_conn->N_adj[face_id];
	double nf[2];
	double lf = nb_partition_edge_get_normal(part, face_id, nf);

	if (N > 1) {
		boundary_mask[face_id] = 0;
		get_internal_face_strain(face_id, face_elems_conn, part,
					 disp, N_qp, lf, nf, strain);
	} else {
		boundary_mask[face_id] = 1;
		get_boundary_face_strain(face_id, face_elems_conn, part,
					 bcond, disp, N_qp, lf, strain);
	}
}

static void get_internal_face_strain(uint32_t face_id,
				     const nb_graph_t *face_elems_conn,
				     const nb_partition_t *const part,
				     const double *disp, uint8_t N_qp,
				     double length, double nf[2],
				     double *strain)
{
	uint16_t N = face_elems_conn->N_adj[face_id];
	uint32_t *adj = face_elems_conn->adj[face_id];
	uint32_t memsize = (3 * N_qp + 2 * N) * sizeof(double);
	char *memblock = NB_SOFT_MALLOC(memsize);
	double *wqp = (void*) memblock;
	double *xqp = (void*) (memblock + N_qp * sizeof(double));
	double *grad_phi = (void*) (memblock + 3 * N_qp * sizeof(double));

	get_quadrature_points(part, face_id, length, N_qp, xqp, wqp);

	for (uint8_t q = 0; q < N_qp; q++) {
		interpolators_eval_grad(part, face_id, nf, N, adj,
					&(xqp[q*2]), grad_phi);
		uint32_t id = face_id * N_qp + q;
		memset(&(strain[id*3]), 0, 3 * sizeof(*strain));
		for (uint16_t i = 0; i < N; i++) {
			uint32_t elem_id = adj[i];
			double u = disp[elem_id * 2];
			double v = disp[elem_id*2+1];
			double dphi_dx = grad_phi[i * 2];
			double dphi_dy = grad_phi[i*2+1];
			strain[id * 3] += dphi_dx * u;
			strain[id*3+1] += dphi_dy * v;
			strain[id*3+2] += dphi_dy * u + dphi_dx * v;
		}
	}
	NB_SOFT_FREE(memsize, memblock);
}

static void get_boundary_face_strain(uint32_t face_id,
				     const nb_graph_t *face_elems_conn,
				     const nb_partition_t *const part,
				     const nb_bcond_t *const bcond,
				     const double *disp, uint8_t N_qp,
				     double length,
				     double *strain)
{
	memset(&(strain[face_id * N_qp * 3]), 0, 3 * N_qp * sizeof(*strain));
}

void nb_cvfa_compute_stress_from_strain(const nb_partition_t *part,
					const nb_material_t *const material,
					nb_analysis2D_t analysis2D,
					uint8_t N_quadrature_points,
					const double* strain,
					double* stress /* Output */)
{
	uint32_t N_faces = nb_partition_get_N_edges(part);
	for (uint32_t i = 0; i < N_faces; i++) {
		double D[4];
		nb_pde_get_constitutive_matrix(D, material, analysis2D);
		
		for (uint8_t q = 0; q < N_quadrature_points; q++) {
			uint32_t id = i * N_quadrature_points + q;
			stress[id * 3] = strain[id * 3] * D[0] +
				strain[id*3+1] * D[1];
			stress[id*3+1] = strain[id * 3] * D[1] +
				strain[id*3+1] * D[2];
			stress[id*3+2] = strain[id*3+2] * D[3];
		}
	}
}
