#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "nb/math_bot.h"
#include "nb/memory_bot.h"
#include "nb/solver_bot.h"
#include "nb/geometric_bot.h"
#include "nb/graph_bot.h"
#include "nb/pde_bot/material.h"
#include "nb/pde_bot/common_solid_mechanics/analysis2D.h"
#include "nb/pde_bot/common_solid_mechanics/formulas.h"
#include "nb/pde_bot/boundary_conditions/bcond.h"
#include "nb/pde_bot/boundary_conditions/bcond_iter.h"

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
		uint32_t id1 = nb_mesh2D_edge_get_1n((prt), (id));	\
		uint32_t id2 = nb_mesh2D_edge_get_2n((prt), (id));	\
									\
		s1[0] = nb_mesh2D_node_get_x((prt), id1);		\
		s1[1] = nb_mesh2D_node_get_y((prt), id1);		\
									\
		s2[0] = nb_mesh2D_node_get_x((prt), id2);		\
		s2[1] = nb_mesh2D_node_get_y((prt), id2);		\
	} while(0)

static uint32_t get_cvfa_memsize(uint32_t N_elems);
static void distribute_cvfa_memory(char *memblock,
				   uint32_t N_elems, double **F,
				   nb_graph_t **face_elems_conn);
static void init_global_matrix(const nb_mesh2D_t *part, nb_sparse_t **K);
static void load_face_elems_conn(const nb_mesh2D_t *part,
				 nb_graph_t *face_elems_conn);
static uint32_t get_N_total_face_adj(const nb_mesh2D_t *part);
static uint16_t face_get_N_ngb(const nb_mesh2D_t *part,
			       uint32_t elem_id1, uint32_t elem_id2);
static uint16_t get_N_ngb_around_right_vtx(const nb_mesh2D_t *part,
					   uint32_t elem_id1,
					   uint32_t elem_id2);
static char *set_elemental_faces(nb_graph_t *face_elems_conn,
				 const nb_mesh2D_t *part,
				 char *memblock, uint32_t elem_id);
static char *set_boundary_face(nb_graph_t *face_elems_conn,
			       const nb_mesh2D_t *part, 
			       char *memblock, uint32_t face_id,
			       uint32_t elem_id);
static char *set_internal_face(nb_graph_t *face_elems_conn,
			       const nb_mesh2D_t *part,
			       char *memblock, uint32_t face_id,
			       uint32_t elem_id1, uint32_t elem_id2);
static void set_local_ngb(uint32_t *adj, const nb_mesh2D_t *part,
			  uint32_t elem_id1, uint32_t elem_id2);
static uint16_t get_ngb_around_right_vtx(const nb_mesh2D_t *part,
					 uint32_t *ngb, uint16_t current_id,
					 uint32_t elem_id1, uint32_t elem_id2);
static void assemble_global_forces(double *F,
				   const nb_mesh2D_t *const part,
				   const nb_material_t *material,
				   bool enable_self_weight,
				   double gravity[2]);
static void integrate_elem_force(const nb_mesh2D_t *part,
				 const nb_material_t *material,
				 bool enable_self_weight,
				 double gravity[2],
				 uint32_t elem_id,
				 double *F);
static void assemble_global_stiffness(nb_sparse_t *K,
				      const nb_mesh2D_t *const part,
				      const nb_graph_t *face_elems_conn,
				      const nb_material_t *material,
				      nb_analysis2D_t analysis2D,
				      nb_analysis2D_params *params2D);
static void assemble_face(uint32_t face_id, nb_sparse_t *K,
			  const nb_mesh2D_t *const part,
			  const nb_graph_t *face_elems_conn,
			  const nb_material_t *material,
			  nb_analysis2D_t analysis2D,
			  nb_analysis2D_params *params2D);
static void face_get_normal(const nb_mesh2D_t *const part,
			    const nb_graph_t *face_elems_conn,
			    uint32_t face_id, double nf[2]);
static uint8_t add_subface_if_intersected(const nb_mesh2D_t *part,
					  uint32_t face_id,
					  const nb_msh3trg_t *msh3,
					  const uint32_t *input_msh3_id,
					  uint32_t trg_id,
					  subface_t *subfaces,
					  uint16_t subface_id);
static void add_face_pairwise(const nb_mesh2D_t *part, uint32_t face_id,
			      subface_t *subfaces, uint16_t subface_id);
static void add_subface_in_closest_trg(const nb_mesh2D_t *part, uint32_t face_id,
				 subface_t *subfaces, uint16_t subface_id);

static void get_face_vtx_outside_msh3(subface_t *subfaces, uint16_t N_sf,
				      const double s1[2], const double s2[2],
				      double alone[2]);
static uint16_t get_face_closest_intersection_to_msh3(subface_t *subfaces,
						      uint16_t N_sf,
						      const double alone[2],
						      double p[2]);
static void assemble_internal_face(uint32_t face_id, nb_sparse_t *K,
				   const nb_mesh2D_t *const part,
				   const nb_graph_t *face_elems_conn,
				   const double D[4], const double nf[2],
				   nb_analysis2D_params *params2D);
static void integrate_Kf(const nb_mesh2D_t *const part, uint32_t face_id,
			 const double D[4], const double nf[2],
			 uint16_t N,  const uint32_t *adj,
			 nb_analysis2D_params *params2D, double *Kf);
static uint16_t load_subfaces(subface_t *subfaces,
			      const nb_mesh2D_t *const part,
			      uint16_t N, const uint32_t *adj,
			      uint32_t face_id);
static void load_msh3trg(const nb_mesh2D_t *const part, uint16_t N,
			 const uint32_t *adj, nb_msh3trg_t *msh3);
static void load_msh3trg_input_vtx_permutation(const nb_msh3trg_t *msh3,
					       uint32_t *id);
static bool face_intersects_trg(const nb_mesh2D_t *part, uint32_t face_id,
				nb_msh3trg_t *msh3, uint32_t trg_id);
static void integrate_subface_in_trg(const nb_mesh2D_t *const part,
				     uint32_t face_id, const double D[4],
				     const double nf[2], uint16_t N,
				     const uint32_t *adj,
				     nb_analysis2D_params *params2D,
				     const subface_t *subface, double *Kf);
static void subface_get_triangle_points(const subface_t *subface,
					const nb_mesh2D_t *part,
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
static void integrate_subface_pairwise(const nb_mesh2D_t *const part,
				       uint32_t face_id, const double D[4],
				       const double nf[2], uint16_t N,
				       const uint32_t *adj,
				       nb_analysis2D_params *params2D,
				       const subface_t *subface, double *Kf);
static void subface_get_grad_pairwise(const double c1[2], const double c2[2],
				      double grad[2]);
static void add_Kf_to_K(uint16_t N, const uint32_t *adj,
			const double *Kf, nb_sparse_t *K);
static int solver(const nb_sparse_t *const A,
		  const double *const b, double* x);
static void get_permutation(const nb_sparse_t *const A,
			    uint32_t *perm, uint32_t *iperm);
static void vector_permutation(uint32_t N, const double *v,
			       const uint32_t *perm, double *vp);
static void compute_strain(double *strain, char *boundary_mask,
			   const nb_graph_t *face_elems_conn,
			   const nb_mesh2D_t *const part,
			   const nb_bcond_t *const bcond,
			   const double *disp);
static void get_face_strain(uint32_t face_id,
			    const nb_graph_t *face_elems_conn,
			    const nb_mesh2D_t *const part,
			    const nb_bcond_t *const bcond,
			    const double *disp,
			    double *strain,
			    char *boundary_mask);
static void get_internal_face_strain(uint32_t face_id,
				     const nb_graph_t *face_elems_conn,
				     const nb_mesh2D_t *const part,
				     const double *disp,
				     double length, double nf[2],
				     double *strain);

static void subface_sum_strain_in_trg(const nb_mesh2D_t *const part,
				      uint32_t face_id, uint16_t N,
				      const uint32_t *adj,
				      const subface_t *subface,
				      const double *disp, double *strain);
static void subface_sum_strain_pairwise(const nb_mesh2D_t *const part,
					uint32_t face_id, uint16_t N,
					const uint32_t *adj,
					const subface_t *subface,
					const double *disp, double *strain);
static void get_boundary_face_strain(uint32_t face_id,
				     const nb_graph_t *face_elems_conn,
				     const nb_mesh2D_t *const part,
				     const nb_bcond_t *const bcond,
				     const double *disp,
				     double length, double *strain);

int nb_cvfa_compute_2D_Solid_Mechanics
			(const nb_mesh2D_t *const part,
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
	uint32_t N_elems = nb_mesh2D_get_N_elems(part);
	uint32_t memsize = get_cvfa_memsize(N_elems);
	char *memblock = nb_soft_allocate_mem(memsize);
	double *F;
	nb_graph_t *face_elems_conn;
	distribute_cvfa_memory(memblock, N_elems, &F, &face_elems_conn);

	nb_sparse_t *K;
	init_global_matrix(part, &K);

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
	nb_sparse_destroy(K);
	nb_graph_finish(face_elems_conn);
	nb_soft_free_mem(memsize, memblock);
	return status;
}

static uint32_t get_cvfa_memsize(uint32_t N_elems)
{
	uint32_t graph_size = nb_graph_get_memsize();
	uint32_t system_size = 2 * N_elems * sizeof(double);
	return graph_size + system_size;
}

static void distribute_cvfa_memory(char *memblock,
				   uint32_t N_elems, double **F,
				   nb_graph_t **face_elems_conn)
{
	uint32_t system_size = 2 * N_elems * sizeof(double);
	*F = (void*) memblock;
	*face_elems_conn = (void*) (memblock + system_size);
}

static void init_global_matrix(const nb_mesh2D_t *part, nb_sparse_t **K)
{
	uint32_t memsize = nb_graph_get_memsize();
	nb_graph_t *graph = nb_soft_allocate_mem(memsize);

	nb_graph_init(graph);
	nb_mesh2D_load_graph(part, graph, NB_ELEMS_LINKED_BY_NODES);

	*K = nb_sparse_create(graph, NULL, 2);

	nb_graph_finish(graph);
	nb_soft_free_mem(memsize, graph);
}

static void load_face_elems_conn(const nb_mesh2D_t *part,
				 nb_graph_t *face_elems_conn)
{
	uint32_t N = nb_mesh2D_get_N_edges(part);
	uint32_t N_total_adj = get_N_total_face_adj(part);
	face_elems_conn->N = N;
	uint32_t memsize = N * (sizeof(*(face_elems_conn->N_adj)) +
				sizeof(*(face_elems_conn->adj))) +
		N_total_adj * sizeof(**(face_elems_conn->adj));
	char *memblock = nb_allocate_mem(memsize);

	face_elems_conn->N_adj = (void*) memblock;
	face_elems_conn->adj = (void*)
		(memblock + N * sizeof(*(face_elems_conn->N_adj)));
	
	uint32_t N_elems = nb_mesh2D_get_N_elems(part);
	memblock += N * (sizeof(*(face_elems_conn->N_adj)) +
			 sizeof(*(face_elems_conn->adj)));

	for (uint32_t i = 0; i < N_elems; i++)		
		memblock = set_elemental_faces(face_elems_conn, part,
					       memblock, i);
}

static uint32_t get_N_total_face_adj(const nb_mesh2D_t *part)
{
	uint32_t N = 0;
	uint32_t N_elems = nb_mesh2D_get_N_elems(part);
	for (uint32_t i = 0; i < N_elems; i++) {
		uint16_t N_adj = nb_mesh2D_elem_get_N_adj(part, i);
		for (uint16_t j = 0; j < N_adj; j++) {
			uint32_t ngb_id = 
				nb_mesh2D_elem_get_ngb(part, i, j);
			if (N_elems <= ngb_id)
				N += 1;
			else if (i < ngb_id)
				N += face_get_N_ngb(part, i, ngb_id);
		}
	}
	return N;
}

static uint16_t face_get_N_ngb(const nb_mesh2D_t *part,
			       uint32_t elem_id1, uint32_t elem_id2)
{
	return 2;
}

static char *set_elemental_faces(nb_graph_t *face_elems_conn,
				 const nb_mesh2D_t *part,
				 char *memblock, uint32_t elem_id)
{
	uint32_t N_elems = nb_mesh2D_get_N_elems(part);
	uint16_t N_adj = nb_mesh2D_elem_get_N_adj(part, elem_id);
	for (uint32_t j = 0; j < N_adj; j++) {
		uint32_t ngb_id = nb_mesh2D_elem_get_ngb(part, elem_id, j);
		if (N_elems <= ngb_id) {
			uint32_t face_id = 
				nb_mesh2D_elem_find_edge(part, elem_id, j);
			memblock = set_boundary_face(face_elems_conn, part,
						     memblock, face_id, elem_id);
		} else if (elem_id < ngb_id) {
			uint32_t face_id =
				nb_mesh2D_elem_find_edge(part, elem_id, j);
			memblock = set_internal_face(face_elems_conn, part,
						     memblock, face_id,
						     elem_id, ngb_id);
		}
	}
	return memblock;
}

static char *set_boundary_face(nb_graph_t *face_elems_conn,
			       const nb_mesh2D_t *part, 
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
			       const nb_mesh2D_t *part,
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

static void set_local_ngb(uint32_t *adj, const nb_mesh2D_t *part,
			  uint32_t elem_id1, uint32_t elem_id2)
{
	adj[0] = elem_id1;
	adj[1] = elem_id2;
}

static void assemble_global_forces(double *F,
				   const nb_mesh2D_t *const part,
				   const nb_material_t *material,
				   bool enable_self_weight,
				   double gravity[2])
{
	uint32_t N_elems = nb_mesh2D_get_N_elems(part);
	memset(F, 0, N_elems * 2 * sizeof(*F));
	for (uint32_t i = 0; i < N_elems; i++) {
		integrate_elem_force(part, material, enable_self_weight,
				     gravity, i, F);
	}
}

static void integrate_elem_force(const nb_mesh2D_t *part,
				 const nb_material_t *material,
				 bool enable_self_weight,
				 double gravity[2],
				 uint32_t elem_id,
				 double *F)
{
	if (enable_self_weight) {
		double area = nb_mesh2D_elem_get_area(part, elem_id);
		double mass = area * nb_material_get_density(material);
		F[elem_id * 2] += mass * gravity[0];
		F[elem_id*2+1] += mass * gravity[1];
	}
}

static void assemble_global_stiffness(nb_sparse_t *K,
				      const nb_mesh2D_t *const part,
				      const nb_graph_t *face_elems_conn,
				      const nb_material_t *material,
				      nb_analysis2D_t analysis2D,
				      nb_analysis2D_params *params2D)
{
	nb_sparse_reset(K);
	uint32_t N_faces = nb_mesh2D_get_N_edges(part);
	for (uint32_t i = 0; i < N_faces; i++) {
		assemble_face(i, K, part, face_elems_conn, material,
			      analysis2D, params2D);
	}
}

static void assemble_face(uint32_t face_id, nb_sparse_t *K,
			  const nb_mesh2D_t *const part,
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

static void face_get_normal(const nb_mesh2D_t *const part,
			    const nb_graph_t *face_elems_conn,
			    uint32_t face_id, double nf[2])
{
	nb_mesh2D_edge_get_normal(part, face_id, nf);

	uint16_t N = face_elems_conn->N_adj[face_id];
	uint32_t N_elems = nb_mesh2D_get_N_elems(part);

	uint32_t elem1 = face_elems_conn->adj[face_id][0];
	uint32_t elem2 = (1 < N)?(face_elems_conn->adj[face_id][1]):(N_elems);

	uint16_t local_fid =
		nb_mesh2D_elem_ngb_get_face(part, elem1, elem2);
	uint32_t v1 = nb_mesh2D_elem_get_adj(part, elem1, local_fid);
	if (v1 != nb_mesh2D_edge_get_1n(part, face_id)) {
		nf[0] *= -1;
		nf[1] *= -1;
	}
}

static void assemble_internal_face(uint32_t face_id, nb_sparse_t *K,
				   const nb_mesh2D_t *const part,
				   const nb_graph_t *face_elems_conn,
				   const double D[4], const double nf[2],
				   nb_analysis2D_params *params2D)
{
	uint16_t N = 2;
	uint32_t *adj = face_elems_conn->adj[face_id];

	double Kf[8];
	integrate_Kf(part, face_id, D, nf, N, adj, params2D, Kf);
	add_Kf_to_K(N, adj, Kf, K);
}

static void integrate_Kf(const nb_mesh2D_t *const part, uint32_t face_id,
			 const double D[4], const double nf[2],
			 uint16_t N,  const uint32_t *adj,
			 nb_analysis2D_params *params2D, double *Kf)
{
	uint32_t sf_size = N * sizeof(subface_t);
	uint32_t memsize = sf_size;
	char *memblock = nb_soft_allocate_mem(memsize);
	subface_t *subfaces = (void*) memblock;

	uint16_t N_sf = load_subfaces(subfaces, part, N, adj, face_id);

	integrate_subface_pairwise(part, face_id, D,
				   nf, N, adj, params2D,
				   subfaces, Kf);	
	nb_soft_free_mem(memsize, memblock);
}

static uint16_t load_subfaces(subface_t *subfaces,
			      const nb_mesh2D_t *const part,
			      uint16_t N, const uint32_t *adj,
			      uint32_t face_id)
{
	add_face_pairwise(part, face_id, subfaces, 0);
	return 1; 
}


static void add_face_pairwise(const nb_mesh2D_t *part, uint32_t face_id,
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

static void integrate_subface_pairwise(const nb_mesh2D_t *const part,
				       uint32_t face_id, const double D[4],
				       const double nf[2], uint16_t N,
				       const uint32_t *adj,
				       nb_analysis2D_params *params2D,
				       const subface_t *subface, double *Kf)
{
	double c1[2], c2[2];
	c1[0] = nb_mesh2D_elem_get_x(part, adj[0]);
	c1[1] = nb_mesh2D_elem_get_y(part, adj[0]);
	c2[0] = nb_mesh2D_elem_get_x(part, adj[1]);
	c2[1] = nb_mesh2D_elem_get_y(part, adj[1]);

	double lf = nb_utils2D_get_dist(subface->xp ,&(subface->xp[2]));
	double factor = lf * params2D->thickness;
	for (uint8_t i = 0; i < 2; i++) {
		double grad[2];
		if (0 == i)
			subface_get_grad_pairwise(c1, c2, grad);
		else
			subface_get_grad_pairwise(c2, c1, grad);
		double Kfi[4];
		subface_get_nodal_contribution(D, nf, grad, Kfi);
		Kf[i * 2] = factor * Kfi[0];
		Kf[i*2+1] = factor * Kfi[1];
		Kf[2 * N + i * 2] = factor * Kfi[2];
		Kf[2 * N + i*2+1] = factor * Kfi[3];
	}
}

static void subface_get_grad_pairwise(const double c1[2], const double c2[2],
				      double grad[2])
{
	double xdiff = c1[0] - c2[0];
	double ydiff = c1[1] - c2[1];
	double dist = nb_utils2D_get_dist(c1, c2);
	double nc[2];
	nc[0] = -xdiff / dist;
	nc[1] = -ydiff / dist;
	double denom = xdiff * nc[0] + ydiff * nc[1];

	grad[0] = nc[0] / denom;
	grad[1] = nc[1] / denom;
}

static void add_Kf_to_K(uint16_t N, const uint32_t *adj,
			const double *Kf, nb_sparse_t *K)
{
	uint16_t size = 2 * N;
	uint32_t i = adj[0];
	uint32_t j = adj[1];
	for (uint8_t m = 0; m < N; m++) {
		uint32_t k = adj[m];
		nb_sparse_add(K, i * 2, k * 2, -Kf[m * 2]);
		nb_sparse_add(K, i * 2, k*2+1, -Kf[m*2+1]);
		nb_sparse_add(K, i*2+1, k * 2, -Kf[size + m * 2]);
		nb_sparse_add(K, i*2+1, k*2+1, -Kf[size + m*2+1]);

		nb_sparse_add(K, j * 2, k * 2, Kf[m * 2]);
		nb_sparse_add(K, j * 2, k*2+1, Kf[m*2+1]);
		nb_sparse_add(K, j*2+1, k * 2, Kf[size + m * 2]);
		nb_sparse_add(K, j*2+1, k*2+1, Kf[size + m*2+1]);
	}
}

static int solver(const nb_sparse_t *const A,
		  const double *const b, double* x)
{
	uint32_t N = nb_sparse_get_size(A);
	uint32_t memsize = 2 * N * (sizeof(uint32_t) + sizeof(double));
	char *memblock = nb_soft_allocate_mem(memsize);
	uint32_t *perm = (void*) memblock;
	uint32_t *iperm = (void*) (memblock + N * sizeof(uint32_t));
	double *br = (void*) (memblock + 2 * N * sizeof(uint32_t));
	double *xr = (void*) (memblock + 2 * N * sizeof(uint32_t) +
			      N * sizeof(double));

	double asym = nb_sparse_get_asym(A);/* TEMPORAL */
	printf("-- K usym: %e\n", asym);     /* TEMPORAL */

	get_permutation(A, perm, iperm);

	nb_sparse_t *Ar = nb_sparse_create_permutation(A, perm, iperm);
	vector_permutation(N, b, perm, br);

	int status = nb_sparse_solve_using_LU(Ar, br, xr, 1);

	vector_permutation(N, xr, iperm, x);
	
	nb_sparse_destroy(Ar);
	nb_soft_free_mem(memsize, memblock);
	return status;
}

static void get_permutation(const nb_sparse_t *const A,
			    uint32_t *perm, uint32_t *iperm)
{
	uint16_t memsize = nb_graph_get_memsize();
	nb_graph_t *graph = nb_soft_allocate_mem(memsize);
	nb_graph_init(graph);
	nb_sparse_get_graph(A, graph);
	nb_graph_labeling(graph, perm, iperm, NB_LABELING_ND);
	nb_graph_finish(graph);

	nb_soft_free_mem(memsize, graph);
}

static void vector_permutation(uint32_t N, const double *v,
			       const uint32_t *perm, double *vp)
{
	for (uint32_t i = 0; i < N; i++)
		vp[i] = v[perm[i]];
}

static void compute_strain(double *strain, char *boundary_mask,
			   const nb_graph_t *face_elems_conn,
			   const nb_mesh2D_t *const part,
			   const nb_bcond_t *const bcond,
			   const double *disp)
{
	;
}

void nb_cvfa_compute_stress_from_strain(const nb_mesh2D_t *part,
					const nb_material_t *const material,
					nb_analysis2D_t analysis2D,
					const double* strain,
					double* stress /* Output */)
{
	uint32_t N_faces = nb_mesh2D_get_N_edges(part);
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
