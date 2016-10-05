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

#include "../integration_mesh.h"
#include "set_bconditions.h"

typedef struct {
	uint8_t N_int;
	double xp[4];
	uint32_t trg[3];
} subface_t;

#define POW2(a) ((a)*(a))

#define LOAD_TRG_INFO(msh, id)						\
	do {								\
		tid1 = nb_msh3trg_elem_get_adj((msh), (id), 0);		\
		tid2 = nb_msh3trg_elem_get_adj((msh), (id), 1);		\
		tid3 = nb_msh3trg_elem_get_adj((msh), (id), 2);		\
									\
		t1[0] = nb_msh3trg_node_get_x((msh), tid1);		\
		t1[1] = nb_msh3trg_node_get_y((msh), tid1);		\
									\
		t2[0] = nb_msh3trg_node_get_x((msh), tid2);		\
		t2[1] = nb_msh3trg_node_get_y((msh), tid2);		\
									\
		t3[0] = nb_msh3trg_node_get_x((msh), tid3);		\
		t3[1] = nb_msh3trg_node_get_y((msh), tid3);		\
	} while(0)

#define LOAD_FACE_INFO(prt, id)						\
	do {								\
		uint32_t id1 = nb_partition_edge_get_1n((prt), (id));	\
		uint32_t id2 = nb_partition_edge_get_2n((prt), (id));	\
									\
		s1[0] = nb_partition_node_get_x((prt), id1);		\
		s1[1] = nb_partition_node_get_y((prt), id1);		\
									\
		s2[0] = nb_partition_node_get_x((prt), id2);		\
		s2[1] = nb_partition_node_get_y((prt), id2);		\
	} while(0)

static uint32_t get_cvfa_memsize(uint32_t N_elems);
static void distribute_cvfa_memory(char *memblock,
				   uint32_t N_elems, double **F,
				   nb_graph_t **face_elems_conn,
				   nb_partition_t **intmsh);
static void init_global_matrix(const nb_partition_t *intmsh, vcn_sparse_t **K);
static void load_face_elems_conn(const nb_partition_t *part,
				 nb_graph_t *face_elems_conn);
static uint32_t get_N_total_face_adj(const nb_partition_t *part);
static uint16_t face_get_N_ngb(const nb_partition_t *part,
			       uint32_t elem_id1, uint32_t elem_id2);
static uint16_t get_N_ngb_around_right_vtx(const nb_partition_t *part,
					   uint32_t elem_id1,
					   uint32_t elem_id2);
static char *set_elemental_faces(nb_graph_t *face_elems_conn,
				 const nb_partition_t *part,
				 char *memblock, uint32_t elem_id);
static char *set_boundary_face(nb_graph_t *face_elems_conn,
			       const nb_partition_t *part, 
			       char *memblock, uint32_t face_id,
			       uint32_t elem_id);
static char *set_internal_face(nb_graph_t *face_elems_conn,
			       const nb_partition_t *part,
			       char *memblock, uint32_t face_id,
			       uint32_t elem_id1, uint32_t elem_id2);
static void set_local_ngb(uint32_t *adj, const nb_partition_t *part,
			  uint32_t elem_id1, uint32_t elem_id2);
static uint16_t get_ngb_around_right_vtx(const nb_partition_t *part,
					 uint32_t *ngb, uint16_t current_id,
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
				      nb_analysis2D_params *params2D);
static void assemble_face(uint32_t face_id, vcn_sparse_t *K,
			  const nb_partition_t *const part,
			  const nb_graph_t *face_elems_conn,
			  const nb_material_t *material,
			  nb_analysis2D_t analysis2D,
			  nb_analysis2D_params *params2D);
static void face_get_normal(const nb_partition_t *const part,
			    const nb_graph_t *face_elems_conn,
			    uint32_t face_id, double nf[2]);
static uint8_t add_subface_if_intersected(const nb_partition_t *part,
					  uint32_t face_id,
					  const nb_msh3trg_t *msh3,
					  const uint32_t *input_msh3_id,
					  uint32_t trg_id,
					  subface_t *subfaces,
					  uint16_t subface_id);
static void add_face_pairwise(const nb_partition_t *part, uint32_t face_id,
			      subface_t *subfaces, uint16_t subface_id);
static void add_subface_in_closest_trg(const nb_partition_t *part, uint32_t face_id,
				 subface_t *subfaces, uint16_t subface_id);

static void get_face_vtx_outside_msh3(subface_t *subfaces, uint16_t N_sf,
				      const double s1[2], const double s2[2],
				      double alone[2]);
static uint16_t get_face_closest_intersection_to_msh3(subface_t *subfaces,
						      uint16_t N_sf,
						      const double alone[2],
						      double p[2]);
static void assemble_internal_face(uint32_t face_id, vcn_sparse_t *K,
				   const nb_partition_t *const part,
				   const nb_graph_t *face_elems_conn,
				   const double D[4], const double nf[2],
				   nb_analysis2D_params *params2D);
static void integrate_Kf(const nb_partition_t *const part, uint32_t face_id,
			 const double D[4], const double nf[2],
			 uint16_t N,  const uint32_t *adj,
			 nb_analysis2D_params *params2D, double *Kf);
static uint16_t load_subfaces(subface_t *subfaces,
			      const nb_partition_t *const part,
			      uint16_t N, const uint32_t *adj,
			      uint32_t face_id);
static void load_msh3trg(const nb_partition_t *const part, uint16_t N,
			 const uint32_t *adj, nb_msh3trg_t *msh3);
static void load_msh3trg_input_vtx_permutation(const nb_msh3trg_t *msh3,
					       uint32_t *id);
static bool face_intersects_trg(const nb_partition_t *part, uint32_t face_id,
				nb_msh3trg_t *msh3, uint32_t trg_id);
static void integrate_subface_in_trg(const nb_partition_t *const part,
				     uint32_t face_id, const double D[4],
				     const double nf[2], uint16_t N,
				     const uint32_t *adj,
				     nb_analysis2D_params *params2D,
				     const subface_t *subface, double *Kf);
static void subface_get_triangle_points(const subface_t *subface,
					const nb_partition_t *part,
					const uint32_t *adj, double t1[2],
					double t2[2], double t3[2]);
static double subface_get_inverse_jacobian(const double t1[2],
					   const double t2[2],
					   const double t3[2],
					   double iJ[4]);
static double subface_get_normalized_length(const subface_t *subface,
					    const double t1[2],
					    const double t2[2],
					    const double t3[2]);
static void get_normalized_point(const double x1[2], const double x2[2],
				 const double x3[2], const double xq[2],
				 double psi[2]);
static void subface_get_normalized_grad(uint8_t i, double grad_xi[2]);
static void subface_get_grad(const double iJ[4], const double grad_xi[2],
			     double grad[2]);
static void subface_get_nodal_contribution(const double D[4],
					   const double nf[2],
					   const double grad[2],
					   double Kfi[4]);
static void integrate_subface_pairwise(const nb_partition_t *const part,
				       uint32_t face_id, const double D[4],
				       const double nf[2], uint16_t N,
				       const uint32_t *adj,
				       nb_analysis2D_params *params2D,
				       const subface_t *subface, double *Kf);
static void subface_get_grad_pairwise(const double c1[2], const double c2[2],
				      double grad[2]);
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
			   const double *disp);
static void get_face_strain(uint32_t face_id,
			    const nb_graph_t *face_elems_conn,
			    const nb_partition_t *const part,
			    const nb_bcond_t *const bcond,
			    const double *disp,
			    double *strain,
			    char *boundary_mask);
static void get_internal_face_strain(uint32_t face_id,
				     const nb_graph_t *face_elems_conn,
				     const nb_partition_t *const part,
				     const double *disp,
				     double length, double nf[2],
				     double *strain);

static void subface_sum_strain_in_trg(const nb_partition_t *const part,
				      uint32_t face_id, uint16_t N,
				      const uint32_t *adj,
				      const subface_t *subface,
				      const double *disp, double *strain);
static void subface_sum_strain_pairwise(const nb_partition_t *const part,
					uint32_t face_id, uint16_t N,
					const uint32_t *adj,
					const subface_t *subface,
					const double *disp, double *strain);
static void get_boundary_face_strain(uint32_t face_id,
				     const nb_graph_t *face_elems_conn,
				     const nb_partition_t *const part,
				     const nb_bcond_t *const bcond,
				     const double *disp,
				     double length, double *strain);

int nb_cvfa_compute_2D_Solid_Mechanics
			(const nb_partition_t *const part,
			 const nb_material_t *const material,
			 const nb_bcond_t *const bcond,
			 bool enable_self_weight, double gravity[2],
			 nb_analysis2D_t analysis2D,
			 nb_analysis2D_params *params2D,
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
	nb_partition_t *intmsh;
	distribute_cvfa_memory(memblock, N_elems, &F, &face_elems_conn,
			       &intmsh);

	nb_cvfa_init_integration_mesh(intmsh);
	nb_cvfa_load_integration_mesh(part, intmsh);
	
	vcn_sparse_t *K;
	init_global_matrix(intmsh, &K);

	nb_graph_init(face_elems_conn);
	load_face_elems_conn(part, face_elems_conn);

	assemble_global_forces(F, part, material, enable_self_weight,
			       gravity);

	assemble_global_stiffness(K, part, face_elems_conn, material,
				  analysis2D, params2D);
	
	nb_cvfa_set_bconditions(part, material, analysis2D, 
				K, F, bcond, 1.0);

	int solver_status = solver(K, F, displacement);
	if (0 != solver_status) {
		status = 1;
		goto CLEANUP_LINEAR_SYSTEM;
	}

	compute_strain(strain, boundary_mask, face_elems_conn, part,
		       bcond, displacement);

	status = 0;
CLEANUP_LINEAR_SYSTEM:
	vcn_sparse_destroy(K);
	nb_graph_finish(face_elems_conn);
	nb_partition_finish(intmsh);
	NB_SOFT_FREE(memsize, memblock);
	return status;
}

static uint32_t get_cvfa_memsize(uint32_t N_elems)
{
	uint32_t graph_size = nb_graph_get_memsize();
	uint32_t system_size = 2 * N_elems * sizeof(double);
	uint32_t intmsh_size = nb_cvfa_get_integration_mesh_memsize();
	return graph_size + system_size + intmsh_size;
}

static void distribute_cvfa_memory(char *memblock,
				   uint32_t N_elems, double **F,
				   nb_graph_t **face_elems_conn,
				   nb_partition_t **intmsh)
{
	uint32_t system_size = 2 * N_elems * sizeof(double);
	uint32_t graph_size = nb_graph_get_memsize();
	*F = (void*) memblock;
	*face_elems_conn = (void*) (memblock + system_size);
	*intmsh = (void*) (memblock + system_size + graph_size);
}

static void init_global_matrix(const nb_partition_t *intmsh, vcn_sparse_t **K)
{
	uint32_t memsize = nb_graph_get_memsize();
	nb_graph_t *graph = NB_SOFT_MALLOC(memsize);

	nb_graph_init(graph);
	nb_partition_load_graph(intmsh, graph, NB_NODES_LINKED_BY_ELEMS);
	nb_graph_extend_adj(graph, 3);

	*K = vcn_sparse_create(graph, NULL, 2);

	nb_graph_finish(graph);
	NB_SOFT_FREE(memsize, graph);
}

static void load_face_elems_conn(const nb_partition_t *part,
				 nb_graph_t *face_elems_conn)
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
	
	uint32_t N_elems = nb_partition_get_N_elems(part);
	memblock += N * (sizeof(*(face_elems_conn->N_adj)) +
			 sizeof(*(face_elems_conn->adj)));

	for (uint32_t i = 0; i < N_elems; i++)		
		memblock = set_elemental_faces(face_elems_conn, part,
					       memblock, i);
}

static uint32_t get_N_total_face_adj(const nb_partition_t *part)
{
	uint32_t N = 0;
	uint32_t N_elems = nb_partition_get_N_elems(part);
	for (uint32_t i = 0; i < N_elems; i++) {
		uint16_t N_adj = nb_partition_elem_get_N_adj(part, i);
		for (uint32_t j = 0; j < N_adj; j++) {
			uint32_t ngb_id = 
				nb_partition_elem_get_ngb(part, i, j);
			if (N_elems <= ngb_id)
				N += 1;
			else if (i < ngb_id)
				N += face_get_N_ngb(part, i, ngb_id);
		}
	}
	return N;
}

static uint16_t face_get_N_ngb(const nb_partition_t *part,
			       uint32_t elem_id1, uint32_t elem_id2)
{
	uint16_t N = 2 +
		get_N_ngb_around_right_vtx(part, elem_id1, elem_id2) +
		get_N_ngb_around_right_vtx(part, elem_id2, elem_id1);
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
				 char *memblock, uint32_t elem_id)
{
	uint32_t N_elems = nb_partition_get_N_elems(part);
	uint16_t N_adj = nb_partition_elem_get_N_adj(part, elem_id);
	for (uint32_t j = 0; j < N_adj; j++) {
		uint32_t ngb_id = nb_partition_elem_get_ngb(part, elem_id, j);
		if (N_elems <= ngb_id) {
			uint32_t face_id = 
				nb_partition_elem_find_edge(part, elem_id, j);
			memblock = set_boundary_face(face_elems_conn, part,
						     memblock, face_id, elem_id);
		} else if (elem_id < ngb_id) {
			uint32_t face_id =
				nb_partition_elem_find_edge(part, elem_id, j);
			memblock = set_internal_face(face_elems_conn, part,
						     memblock, face_id,
						     elem_id, ngb_id);
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
			       char *memblock, uint32_t face_id,
			       uint32_t elem_id1, uint32_t elem_id2)
{
       	uint16_t N_ngb = face_get_N_ngb(part, elem_id1, elem_id2);
	face_elems_conn->N_adj[face_id] = N_ngb;
	face_elems_conn->adj[face_id] = (void*) memblock;
	memblock += N_ngb * sizeof(**(face_elems_conn->adj));

	set_local_ngb(face_elems_conn->adj[face_id],
		      part, elem_id1, elem_id2);

	return memblock;

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
				      nb_analysis2D_params *params2D)
{
	vcn_sparse_reset(K);
	uint32_t N_faces = nb_partition_get_N_edges(part);
	for (uint32_t i = 0; i < N_faces; i++) {
		assemble_face(i, K, part, face_elems_conn, material,
			      analysis2D, params2D);
	}
}

static void assemble_face(uint32_t face_id, vcn_sparse_t *K,
			  const nb_partition_t *const part,
			  const nb_graph_t *face_elems_conn,
			  const nb_material_t *material,
			  nb_analysis2D_t analysis2D,
			  nb_analysis2D_params *params2D)
{
	double D[4];
	nb_pde_get_constitutive_matrix(D, material, analysis2D);

	double nf[2];
	face_get_normal(part, face_elems_conn, face_id, nf);

	uint16_t N = face_elems_conn->N_adj[face_id];
	if (1 < N)
		assemble_internal_face(face_id, K, part, face_elems_conn,
				       D, nf, params2D);
}

static void face_get_normal(const nb_partition_t *const part,
			    const nb_graph_t *face_elems_conn,
			    uint32_t face_id, double nf[2])
{
	nb_partition_edge_get_normal(part, face_id, nf);

	uint16_t N = face_elems_conn->N_adj[face_id];
	uint32_t N_elems = nb_partition_get_N_elems(part);

	uint32_t elem1 = face_elems_conn->adj[face_id][0];
	uint32_t elem2 = (1 < N)?(face_elems_conn->adj[face_id][1]):(N_elems);

	uint16_t local_fid =
		nb_partition_elem_ngb_get_face(part, elem1, elem2);
	uint32_t v1 = nb_partition_elem_get_adj(part, elem1, local_fid);
	if (v1 != nb_partition_edge_get_1n(part, face_id)) {
		nf[0] *= -1;
		nf[1] *= -1;
	}
}

static void assemble_internal_face(uint32_t face_id, vcn_sparse_t *K,
				   const nb_partition_t *const part,
				   const nb_graph_t *face_elems_conn,
				   const double D[4], const double nf[2],
				   nb_analysis2D_params *params2D)
{
	uint16_t N = face_elems_conn->N_adj[face_id];
	uint32_t *adj = face_elems_conn->adj[face_id];

	uint32_t memsize = 4 * N * sizeof(double);
	char* memblock = NB_SOFT_MALLOC(memsize);
	double *Kf = (void*) memblock;

	integrate_Kf(part, face_id, D, nf, N, adj, params2D, Kf);

	add_Kf_to_K(N, adj, Kf, K);

	NB_SOFT_FREE(memsize, memblock);
}

static void integrate_Kf(const nb_partition_t *const part, uint32_t face_id,
			 const double D[4], const double nf[2],
			 uint16_t N,  const uint32_t *adj,
			 nb_analysis2D_params *params2D, double *Kf)
{
	uint32_t sf_size = N * sizeof(subface_t);
	uint32_t memsize = sf_size;
	char *memblock = NB_SOFT_MALLOC(memsize);
	subface_t *subfaces = (void*) memblock;

	uint16_t N_sf = load_subfaces(subfaces, part, N, adj, face_id);
	
	memset(Kf, 0, 4 * N * sizeof(double));
	for (uint16_t i = 0; i < N_sf; i++) {
		if (0 < subfaces[i].N_int)
			integrate_subface_in_trg(part, face_id, D, nf, N, adj,
						 params2D, &(subfaces[i]), Kf);
		else
			integrate_subface_pairwise(part, face_id, D,
						   nf, N, adj, params2D,
						   &(subfaces[i]), Kf);
	}
	
	NB_SOFT_FREE(memsize, memblock);
}

static uint16_t load_subfaces(subface_t *subfaces,
			      const nb_partition_t *const part,
			      uint16_t N, const uint32_t *adj,
			      uint32_t face_id)
{
	uint32_t msh3_size = nb_msh3trg_get_memsize();
	uint32_t memsize = msh3_size + N * sizeof(uint32_t);
	char *memblock = NB_SOFT_MALLOC(memsize);
	nb_msh3trg_t *msh3 = (void*) memblock;
	uint32_t *id = (void*) (memblock + msh3_size);

	nb_msh3trg_init(msh3);
	load_msh3trg(part, N, adj, msh3);
	load_msh3trg_input_vtx_permutation(msh3, id);
	
	uint16_t N_trg = nb_msh3trg_get_N_elems(msh3);
	uint16_t N_sf = 0;
	uint8_t end_trg = 0;
	for (uint16_t i = 0; i < N_trg; i++) {
		uint8_t N_int = add_subface_if_intersected(part, face_id,
							   msh3, id, i,
							   subfaces, N_sf);
		if (0 < N_int)
			N_sf += 1;
		
		if (1 == N_int)
			end_trg += 1;
			
	}

	if (end_trg < 1) {
		add_face_pairwise(part, face_id, subfaces, N_sf);
		N_sf += 1;
	} else if (end_trg < 2) {
		add_subface_in_closest_trg(part, face_id, subfaces, N_sf);
		N_sf += 1;
	}

	nb_msh3trg_finish(msh3);
	NB_SOFT_FREE(memsize, memblock);
	return N_sf;
}

static void load_msh3trg(const nb_partition_t *const part, uint16_t N,
			 const uint32_t *adj, nb_msh3trg_t *msh3)
{
	uint32_t mesh_size = nb_mesh_get_memsize();
	uint16_t vtx_size = 2 * N * sizeof(double);
	uint32_t memsize = mesh_size + vtx_size;
	char *memblock = NB_SOFT_MALLOC(memsize);
	nb_mesh_t *mesh = (void*) memblock;
	double *vtx = (void*) (memblock + mesh_size);
	
	for (uint16_t i = 0; i < N; i++) {
		vtx[i * 2] = nb_partition_elem_get_x(part, adj[i]);
		vtx[i*2+1] = nb_partition_elem_get_y(part, adj[i]);
	}

	nb_mesh_init(mesh);
	nb_mesh_get_delaunay(mesh, N, vtx);
	nb_msh3trg_load_from_mesh(msh3, mesh);
	nb_mesh_finish(mesh);

	NB_SOFT_FREE(memsize, memblock);
}

static void load_msh3trg_input_vtx_permutation(const nb_msh3trg_t *msh3,
					       uint32_t *id)
{
	uint16_t N = nb_msh3trg_get_N_invtx(msh3);
	for (uint16_t i = 0; i < N; i++) {
		uint32_t nid = nb_msh3trg_get_invtx(msh3, i);
		id[nid] = i;
	}
}

static uint8_t add_subface_if_intersected(const nb_partition_t *part,
					  uint32_t face_id,
					  const nb_msh3trg_t *msh3,
					  const uint32_t *input_msh3_id,
					  uint32_t trg_id,
					  subface_t *subfaces,
					  uint16_t subface_id)
{
	uint32_t tid1, tid2, tid3;
	double t1[2], t2[2], t3[2];
	LOAD_TRG_INFO(msh3, trg_id);

	double s1[2], s2[2];
	LOAD_FACE_INFO(part, face_id);

	uint8_t N_int = 0;
	double xp[4], p[2];

	if (vcn_utils2D_are_sgm_intersected(s1, s2, t1, t2, p)) {
		xp[N_int * 2] = p[0];
		xp[N_int*2+1] = p[1];
		N_int += 1;
	}

	if (vcn_utils2D_are_sgm_intersected(s1, s2, t2, t3, p)) {
		xp[N_int * 2] = p[0];
		xp[N_int*2+1] = p[1];
		N_int += 1;
	}
	
	if (2 > N_int) {
		if (vcn_utils2D_are_sgm_intersected(s1, s2, t3, t1, p)) {
			xp[N_int * 2] = p[0];
			xp[N_int*2+1] = p[1];
			N_int += 1;
		}
	}

	if (1 == N_int) {
		if (vcn_utils2D_pnt_lies_in_trg(t1, t2, t3, s1)) {
			xp[2] = s1[0];
			xp[3] = s1[1];
		} else {
			xp[2] = s2[0];
			xp[3] = s2[1];			
		}
	}
	
	if (0 < N_int) {
		subfaces[subface_id].N_int = N_int;
		memcpy(subfaces[subface_id].xp, xp, 4 * sizeof(*xp));
		subfaces[subface_id].trg[0] = input_msh3_id[tid1];
		subfaces[subface_id].trg[1] = input_msh3_id[tid2];
		subfaces[subface_id].trg[2] = input_msh3_id[tid3];
	}
	
	return N_int;
}

static void add_face_pairwise(const nb_partition_t *part, uint32_t face_id,
			      subface_t *subfaces, uint16_t subface_id)
{
	double s1[2], s2[2];
	LOAD_FACE_INFO(part, face_id);

	subfaces[subface_id].N_int = 0;
	memcpy(subfaces[subface_id].xp, s1, 2 * sizeof(*s1));
	memcpy(&(subfaces[subface_id].xp[2]), s2, 2 * sizeof(*s2));
	subfaces[subface_id].trg[0] = 0;
	subfaces[subface_id].trg[1] = 1;
}

static void add_subface_in_closest_trg(const nb_partition_t *part, uint32_t face_id,
				 subface_t *subfaces, uint16_t subface_id)
{
	double s1[2], s2[2];
	LOAD_FACE_INFO(part, face_id);

	double alone[2];
	get_face_vtx_outside_msh3(subfaces, subface_id, s1, s2, alone);

	double p[2];
	uint16_t closest_id =
		get_face_closest_intersection_to_msh3(subfaces, subface_id,
						      alone, p);
	subfaces[subface_id].N_int = 1;
	memcpy(subfaces[subface_id].xp, p, 2 * sizeof(*p));
	memcpy(&(subfaces[subface_id].xp[2]), alone, 2 * sizeof(*alone));
	subfaces[subface_id].trg[0] = subfaces[closest_id].trg[0];
	subfaces[subface_id].trg[1] = subfaces[closest_id].trg[1];
	subfaces[subface_id].trg[2] = subfaces[closest_id].trg[2];
}

static void get_face_vtx_outside_msh3(subface_t *subfaces, uint16_t N_sf,
				      const double s1[2], const double s2[2],
				      double alone[2])
{
	bool s1_is_outside = true;
	for (uint16_t i = 0; i < N_sf; i++) {
		if (1 == subfaces[i].N_int) {
			double e1 = fabs(s1[0] - subfaces[i].xp[2]);
			double e2 = fabs(s1[1] - subfaces[i].xp[3]);
			if (e1 < 1e-20 && e2 < 1e-20)
				s1_is_outside = false;
			break;
		}
	}
	if (s1_is_outside)
		memcpy(alone, s1, 2 * sizeof(*s1));
	else
		memcpy(alone, s2, 2 * sizeof(*s2));
}

static uint16_t get_face_closest_intersection_to_msh3(subface_t *subfaces,
						      uint16_t N_sf,
						      const double alone[2],
						      double p[2])
/* PENDING: SLOW FUNCTION (CALCULATE DIST AND DUPLICATE CHECKS) */
{
	uint16_t sf_id = N_sf;
	double min = 1e30;
	for (uint16_t i = 0; i < N_sf; i++) {
		double d = vcn_utils2D_get_dist(alone, subfaces[i].xp);
		if (d < min) {
			min = d;
			memcpy(p, subfaces[i].xp, 2 * sizeof(*p));
			sf_id = i;
		}
		d = vcn_utils2D_get_dist(alone, &(subfaces[i].xp[2]));
		if (d < min) {
			min = d;
			memcpy(p, &(subfaces[i].xp[2]), 2 * sizeof(*p));
			sf_id = i;
		}
	}
	return sf_id;
}


static void integrate_subface_in_trg(const nb_partition_t *const part,
				     uint32_t face_id, const double D[4],
				     const double nf[2], uint16_t N,
				     const uint32_t *adj,
				     nb_analysis2D_params *params2D,
				     const subface_t *subface, double *Kf)
{
	double t1[2], t2[2], t3[2];
	subface_get_triangle_points(subface, part, adj, t1, t2, t3);

	double iJ[4];
	double detJ = subface_get_inverse_jacobian(t1, t2, t3, iJ);

	double lfn = subface_get_normalized_length(subface, t1, t2, t3);

	double factor = lfn * detJ * params2D->thickness;
	for (uint8_t i = 0; i < 3; i++) {
		double grad_xi[2];
		subface_get_normalized_grad(i, grad_xi);
		double grad[2];
		subface_get_grad(iJ, grad_xi, grad);
		double Kfi[4];
		subface_get_nodal_contribution(D, nf, grad, Kfi);
		uint16_t id = subface->trg[i];
		Kf[id * 2] += factor * Kfi[0];
		Kf[id*2+1] += factor * Kfi[1];
		Kf[2 * N + id * 2] += factor * Kfi[2];
		Kf[2 * N + id*2+1] += factor * Kfi[3];
	}
}

static void subface_get_triangle_points(const subface_t *subface,
					const nb_partition_t *part,
					const uint32_t *adj, double t1[2],
					double t2[2], double t3[2])
{
	uint32_t id1 = adj[subface->trg[0]];
	uint32_t id2 = adj[subface->trg[1]];
	uint32_t id3 = adj[subface->trg[2]];

	t1[0] = nb_partition_elem_get_x(part, id1);
	t1[1] = nb_partition_elem_get_y(part, id1);

	t2[0] = nb_partition_elem_get_x(part, id2);
	t2[1] = nb_partition_elem_get_y(part, id2);

	t3[0] = nb_partition_elem_get_x(part, id3);
	t3[1] = nb_partition_elem_get_y(part, id3);
}

static double subface_get_inverse_jacobian(const double t1[2],
					   const double t2[2],
					   const double t3[2],
					   double iJ[4])
{
	/* Jacobian = D_{psi} x */
	iJ[0] = t2[0] - t1[0];
	iJ[1] = t2[1] - t1[1];
	iJ[2] = t3[0] - t1[0];
	iJ[3] = t3[1] - t1[1];

	double det = vcn_matrix_2X2_inverse_destructive(iJ);

	return det;
}

static double subface_get_normalized_length(const subface_t *subface,
					    const double t1[2],
					    const double t2[2],
					    const double t3[2])
{
	double psi1[2];
	get_normalized_point(t1, t2, t3, subface->xp, psi1);

	double psi2[2];
	get_normalized_point(t1, t2, t3, &(subface->xp[2]), psi2);

	return vcn_utils2D_get_dist(psi1, psi2);
}

static void get_normalized_point(const double x1[2], const double x2[2],
				 const double x3[2], const double xq[2],
				 double psi[2])
{
	double A[4];
	A[0] = x2[0] - x1[0];
	A[1] = x3[0] - x1[0];
	A[2] = x2[1] - x1[1];
	A[3] = x3[1] - x1[1];
	
	double b[2];
	b[0] = xq[0] - x1[0];
	b[1] = xq[1] - x1[1];

	vcn_matrix_2X2_inverse_destructive(A);

	psi[0] = A[0] * b[0] + A[1] * b[1];
	psi[1] = A[2] * b[0] + A[3] * b[1];
}

static void subface_get_normalized_grad(uint8_t i, double grad_xi[2])
{
	if (0 == i) {
		grad_xi[0] = -1;
		grad_xi[1] = -1;
	} else if (1 == i) {
		grad_xi[0] = 1;
		grad_xi[1] = 0;
	} else {
		grad_xi[0] = 0;
		grad_xi[1] = 1;
	}
}

static void subface_get_grad(const double iJ[4], const double grad_xi[2],
			     double grad[2])
{
	grad[0] = iJ[0] * grad_xi[0] + iJ[1] * grad_xi[1];
	grad[1] = iJ[2] * grad_xi[0] + iJ[3] * grad_xi[1];
}

static void subface_get_nodal_contribution(const double D[4],
					   const double nf[2],
					   const double grad[2],
					   double Kfi[4])
{
	Kfi[0] = grad[0] * nf[0] * D[0] + grad[1] * nf[1] * D[3];
	Kfi[1] = grad[1] * nf[0] * D[1] + grad[0] * nf[1] * D[3];
	Kfi[2] = grad[0] * nf[1] * D[1] + grad[1] * nf[0] * D[3];
	Kfi[3] = grad[1] * nf[1] * D[2] + grad[0] * nf[0] * D[3];
}

static void integrate_subface_pairwise(const nb_partition_t *const part,
				       uint32_t face_id, const double D[4],
				       const double nf[2], uint16_t N,
				       const uint32_t *adj,
				       nb_analysis2D_params *params2D,
				       const subface_t *subface, double *Kf)
{
	double c1[2], c2[2];
	c1[0] = nb_partition_elem_get_x(part, adj[0]);
	c1[1] = nb_partition_elem_get_y(part, adj[0]);
	c2[0] = nb_partition_elem_get_x(part, adj[1]);
	c2[1] = nb_partition_elem_get_y(part, adj[1]);

	double lf = vcn_utils2D_get_dist(subface->xp ,&(subface->xp[2]));
	double factor = lf * params2D->thickness;
	for (uint8_t i = 0; i < 2; i++) {
		double grad[2];
		if (0 == i)
			subface_get_grad_pairwise(c1, c2, grad);
		else
			subface_get_grad_pairwise(c2, c1, grad);
		double Kfi[4];
		subface_get_nodal_contribution(D, nf, grad, Kfi);
		Kf[i * 2] += factor * Kfi[0];
		Kf[i*2+1] += factor * Kfi[1];
		Kf[2 * N + i * 2] += factor * Kfi[2];
		Kf[2 * N + i*2+1] += factor * Kfi[3];
	}
}

static void subface_get_grad_pairwise(const double c1[2], const double c2[2],
				      double grad[2])
{
	double xdiff = c1[0] - c2[0];
	double ydiff = c1[1] - c2[1];
	if (fabs(xdiff) < 1e-16) {
		grad[0] = 0;
		grad[1] = 1 / ydiff;
	} else if (fabs(ydiff) < 1e-16) {
		grad[0] = 1 / xdiff;
		grad[1] = 0;
	} else {
		grad[0] = 0.5 / xdiff;
		grad[1] = 0.5 / ydiff;
	}
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
			   const double *disp)
{
	uint32_t N_faces = nb_partition_get_N_edges(part);

 	for (uint32_t i = 0; i < N_faces; i++) {
		get_face_strain(i, face_elems_conn, part, bcond,
				disp, strain, boundary_mask);
	}
}

static void get_face_strain(uint32_t face_id,
			    const nb_graph_t *face_elems_conn,
			    const nb_partition_t *const part,
			    const nb_bcond_t *const bcond,
			    const double *disp,
			    double *strain,
			    char *boundary_mask)
{
	uint16_t N = face_elems_conn->N_adj[face_id];
	double nf[2];
	double lf = nb_partition_edge_get_normal(part, face_id, nf);

	if (N > 1) {
		boundary_mask[face_id] = 0;
		get_internal_face_strain(face_id, face_elems_conn, part,
					 disp, lf, nf, strain);
	} else {
		boundary_mask[face_id] = 1;
		get_boundary_face_strain(face_id, face_elems_conn, part,
					 bcond, disp, lf, strain);
	}
}

static void get_internal_face_strain(uint32_t face_id,
				     const nb_graph_t *face_elems_conn,
				     const nb_partition_t *const part,
				     const double *disp,
				     double length, double nf[2],
				     double *strain)
{
 	uint32_t N = face_elems_conn->N_adj[face_id];
 	uint32_t *adj = face_elems_conn->adj[face_id];
	uint32_t sf_size = N * sizeof(subface_t);
	uint32_t memsize = sf_size;
	char *memblock = NB_SOFT_MALLOC(memsize);
	subface_t *subfaces = (void*) memblock;

	uint16_t N_sf = load_subfaces(subfaces, part, N, adj, face_id);
	
	memset(&(strain[face_id*3]), 0, 3 * sizeof(*strain));
	for (uint16_t i = 0; i < N_sf; i++) {
		if (0 < subfaces[i].N_int)
			subface_sum_strain_in_trg(part, face_id, N, adj,
						  &(subfaces[i]),
						  disp, strain);
		else
			subface_sum_strain_pairwise(part, face_id, N, adj,
						    &(subfaces[i]),
						    disp, strain);
	}
	strain[face_id * 3] /= length;
	strain[face_id*3+1] /= length;
	strain[face_id*3+2] /= length;

	NB_SOFT_FREE(memsize, memblock);
}

static void subface_sum_strain_in_trg(const nb_partition_t *const part,
				      uint32_t face_id, uint16_t N,
				      const uint32_t *adj,
				      const subface_t *subface,
				      const double *disp, double *strain)
{
	double t1[2], t2[2], t3[2];
	subface_get_triangle_points(subface, part, adj, t1, t2, t3);

	double iJ[4];
	double detJ = subface_get_inverse_jacobian(t1, t2, t3, iJ);

	double lfn = subface_get_normalized_length(subface, t1, t2, t3);

	double factor = lfn * detJ;
	for (uint8_t i = 0; i < 3; i++) {
		double grad_xi[2];
		subface_get_normalized_grad(i, grad_xi);
		double grad[2];
		subface_get_grad(iJ, grad_xi, grad);
		uint32_t id = adj[subface->trg[i]];
		strain[face_id * 3] += factor * (grad[0] * disp[id * 2]);
		strain[face_id*3+1] += factor * (grad[1] * disp[id*2+1]);
		strain[face_id*3+2] += factor * (grad[1] * disp[id * 2] +
						 grad[0] * disp[id*2+1]);
	}
}

static void subface_sum_strain_pairwise(const nb_partition_t *const part,
					uint32_t face_id, uint16_t N,
					const uint32_t *adj,
					const subface_t *subface,
					const double *disp, double *strain)
{
	double lf = vcn_utils2D_get_dist(subface->xp ,&(subface->xp[2]));
	
	double c1[2], c2[2];
	c1[0] = nb_partition_elem_get_x(part, adj[0]);
	c1[1] = nb_partition_elem_get_y(part, adj[0]);
	c2[0] = nb_partition_elem_get_x(part, adj[1]);
	c2[1] = nb_partition_elem_get_y(part, adj[1]);

	for (uint8_t i = 0; i < 2; i++) {
		double grad[2];
		if (0 == i)
			subface_get_grad_pairwise(c1, c2, grad);
		else
			subface_get_grad_pairwise(c2, c1, grad);
		strain[face_id * 3] += lf * (grad[0] * disp[adj[i] * 2]);
		strain[face_id*3+1] += lf * (grad[1] * disp[adj[i]*2+1]);
		strain[face_id*3+2] += lf * (grad[1] * disp[adj[i] * 2] +
					     grad[0] * disp[adj[i]*2+1]);
	}
}

static void get_boundary_face_strain(uint32_t face_id,
				     const nb_graph_t *face_elems_conn,
				     const nb_partition_t *const part,
				     const nb_bcond_t *const bcond,
				     const double *disp,
				     double length,
				     double *strain)
{
	memset(&(strain[face_id * 3]), 0, 3 * sizeof(*strain));
}

void nb_cvfa_compute_stress_from_strain(const nb_partition_t *part,
					const nb_material_t *const material,
					nb_analysis2D_t analysis2D,
					const double* strain,
					double* stress /* Output */)
{
	uint32_t N_faces = nb_partition_get_N_edges(part);
	for (uint32_t i = 0; i < N_faces; i++) {
		double D[4];
		nb_pde_get_constitutive_matrix(D, material, analysis2D);
		
		stress[i * 3] = (strain[i * 3] * D[0] +
				 strain[i*3+1] * D[1]);
		stress[i*3+1] = (strain[i * 3] * D[1] +
				 strain[i*3+1] * D[2]);
		stress[i*3+2] = strain[i*3+2] * D[3];
	}
}
