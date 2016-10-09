#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "nb/math_bot.h"
#include "nb/memory_bot.h"
#include "nb/container_bot.h"
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

#define POW2(a) ((a)*(a))

typedef struct subface_s subface_t;

typedef struct {
	double nf[2];
	uint32_t elems[2];
	double x1[2], x2[2];
	uint16_t N_sf;
	subface_t **subfaces;
} face_t;

struct subface_s {
	uint8_t N_int;
	double x1[2], x2[2];
	uint32_t trg_id;
};

static uint32_t get_cvfa_memsize(uint32_t N_elems, uint32_t N_faces);
static void distribute_cvfa_memory(char *memblock, uint32_t N_elems,
				   uint32_t N_faces, double **F,
				   nb_partition_t **intmsh,
				   nb_graph_t **trg_x_vol,
				   face_t ***faces);
static void init_global_matrix(vcn_sparse_t **K, const nb_graph_t *trg_x_vol,
			       const nb_partition_t *intmsh);
static void load_faces(const nb_partition_t *part,
		       const nb_partition_t *intmsh,
		       const nb_graph_t *trg_x_vol,
		       face_t **faces);
static void get_face_elems(const nb_partition_t *part, face_t **faces);
static void define_face_elems(const nb_partition_t *part,
			      face_t **faces, uint32_t elem_id,
			      uint32_t local_face_id);
static void load_subfaces(face_t **faces, uint32_t face_id,
			  const nb_partition_t *const intmsh,
			  const nb_graph_t *trg_x_vol);
static uint8_t add_subface_if_intersected(nb_membank_t *membank,
					  const nb_partition_t *intmsh,
					  const uint32_t *trg_adj,
					  face_t **faces, uint32_t elem_trg_id,
					  uint32_t face_id,
					  nb_container_t *subfaces);
static void add_subface_in_closest_trg(nb_membank_t *membank,
				       const nb_partition_t *intmsh,
				       face_t **faces, uint32_t face_id,
				       nb_container_t *subfaces);

static void get_face_vtx_outside_intmsh(const nb_container_t *subfaces,
					const face_t *face,
					double alone[2]);
static uint16_t get_face_closest_intersection_to_intmsh
					(const nb_container_t *subfaces,
					 const double alone[2],
					 double p[2]);
static void set_subfaces(nb_membank_t *membank, face_t *face,
			 nb_container_t *subfaces);
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
				      const nb_partition_t *intmsh,
				      face_t **faces,
				      const nb_material_t *material,
				      nb_analysis2D_t analysis2D,
				      nb_analysis2D_params *params2D);
static void assemble_face(vcn_sparse_t *K,
			  const nb_partition_t *const part,
			  const nb_partition_t *intmsh,
			  face_t *face,
			  const nb_material_t *material,
			  nb_analysis2D_t analysis2D,
			  nb_analysis2D_params *params2D);
static void assemble_internal_face(vcn_sparse_t *K,
				   const nb_partition_t *const part,
				   const nb_partition_t *intmsh,
				   face_t *face, const double D[4],
				   nb_analysis2D_params *params2D);
static void integrate_subfaces(vcn_sparse_t *K,
			       const nb_partition_t *const part,
			       const nb_partition_t *intmsh,
			       face_t *face, const double D[4],
			       nb_analysis2D_params *params2D);
static void integrate_Kf(const nb_partition_t *const part,
			 const nb_partition_t *intmsh, face_t *face,
			 uint16_t subface_id, const double D[4],
			 nb_analysis2D_params *params2D, double *Kf);
static void load_triangle_points(const nb_partition_t *intmsh,
				 uint32_t trg_id, double t1[2],
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
static void add_Kf_to_K(face_t *face, const nb_partition_t *intmsh,
			uint16_t subface_id, const double *Kf,
			vcn_sparse_t *K);
static void integrate_pairwise(vcn_sparse_t *K,
			       const nb_partition_t *const part,
			       face_t *faces, const double D[4],
			       nb_analysis2D_params *params2D);
static void integrate_Kf_pairwise(const nb_partition_t *const part,
				  face_t *face, const double D[4],
				  nb_analysis2D_params *params2D, double *Kf);
static void face_get_grad_pairwise(const double c1[2], const double c2[2],
				      double grad[2]);
static void add_Kf_to_K_pairwise(face_t *face, const double *Kf,
				 vcn_sparse_t *K);
static int solver(const vcn_sparse_t *const A,
		  const double *const b, double* x);
static void get_permutation(const vcn_sparse_t *const A,
			    uint32_t *perm, uint32_t *iperm);
static void vector_permutation(uint32_t N, const double *v,
			       const uint32_t *perm, double *vp);
static void compute_strain(double *strain, char *boundary_mask,
			   face_t **faces,
			   const nb_partition_t *const part,
			   const nb_partition_t *intmsh,
			   const nb_bcond_t *const bcond,
			   const double *disp);
static void get_face_strain(face_t **faces, uint32_t face_id,
			    const nb_partition_t *const part,
			    const nb_partition_t *intmsh,
			    const nb_bcond_t *const bcond,
			    const double *disp,
			    double *strain,
			    char *boundary_mask);
static void get_internal_face_strain(face_t **faces, uint32_t face_id,
				     const nb_partition_t *const part,
				     const nb_partition_t *intmsh,
				     const double *disp, double *strain);

static void get_subfaces_strain(face_t **faces, uint32_t face_id,
				const nb_partition_t *const part,
				const nb_partition_t *intmsh,
				const double *disp, double *strain);
static void subface_sum_strain_in_trg(const nb_partition_t *const part,
				      const nb_partition_t *intmsh,
				      uint32_t face_id,
				      const subface_t *subface,
				      const double *disp, double *strain);
static void subface_sum_strain_pairwise(const nb_partition_t *const part,
					const nb_partition_t *intmsh,
					uint32_t face_id, 
					const subface_t *subface,
					const double *disp, double *strain);
static void get_pairwise_strain(face_t **faces, uint32_t face_id,
				const nb_partition_t *part,
				const double *disp, double *strain);
static void get_boundary_face_strain(face_t **faces, uint32_t face_id,
				     const nb_partition_t *const part,
				     const nb_bcond_t *bcond,
				     const double *disp, double *strain);
static void finish_faces(uint32_t N_faces, face_t **faces);

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
	uint32_t N_faces = nb_partition_get_N_edges(part);
	uint32_t memsize = get_cvfa_memsize(N_elems, N_faces);
	char *memblock = NB_SOFT_MALLOC(memsize);
	double *F;
	nb_partition_t *intmsh;
	nb_graph_t *trg_x_vol;
	face_t **faces;
	distribute_cvfa_memory(memblock, N_elems, N_faces, &F,
			       &intmsh, &trg_x_vol, &faces);

	nb_cvfa_init_integration_mesh(intmsh);
	nb_cvfa_load_integration_mesh(part, intmsh);

	nb_graph_init(trg_x_vol);
	nb_cvfa_correlate_partition_and_integration_mesh(part, intmsh,
							 trg_x_vol);

	vcn_sparse_t *K;
	init_global_matrix(&K, trg_x_vol, intmsh);

	load_faces(part, intmsh, trg_x_vol, faces);

	assemble_global_forces(F, part, material, enable_self_weight,
			       gravity);

	assemble_global_stiffness(K, part, intmsh, faces, material,
				  analysis2D, params2D);
	
	nb_cvfa_set_bconditions(part, material, analysis2D, 
				K, F, bcond, 1.0);

	int solver_status = solver(K, F, displacement);
	if (0 != solver_status) {
		status = 1;
		goto CLEANUP_LINEAR_SYSTEM;
	}

	compute_strain(strain, boundary_mask, faces, part,
		       intmsh, bcond, displacement);

	status = 0;
CLEANUP_LINEAR_SYSTEM:
	finish_faces(N_faces, faces);
	vcn_sparse_destroy(K);
	nb_graph_finish(trg_x_vol);
	nb_partition_finish(intmsh);
	NB_SOFT_FREE(memsize, memblock);
	return status;
}

static uint32_t get_cvfa_memsize(uint32_t N_elems, uint32_t N_faces)
{
	uint32_t system_size = 2 * N_elems * sizeof(double);
	uint32_t intmsh_size = nb_cvfa_get_integration_mesh_memsize();
	uint32_t graph_size = nb_graph_get_memsize();
	uint32_t faces_size = N_faces * (sizeof(void*) + sizeof(face_t));
	return graph_size + system_size + intmsh_size + faces_size;
}

static void distribute_cvfa_memory(char *memblock, uint32_t N_elems,
				   uint32_t N_faces, double **F,
				   nb_partition_t **intmsh,
				   nb_graph_t **trg_x_vol,
				   face_t ***faces)
{
	uint32_t system_size = 2 * N_elems * sizeof(double);
	uint32_t intmsh_size = nb_cvfa_get_integration_mesh_memsize();
	uint32_t graph_size = nb_graph_get_memsize();
	*F = (void*) memblock;
	*intmsh = (void*) (memblock + system_size);
	*trg_x_vol = (void*) (memblock + system_size + intmsh_size);
	*faces = (void*) (memblock + system_size + intmsh_size +
				    graph_size);
	memblock +=  system_size + intmsh_size + graph_size +
		N_faces * sizeof(void*);
	for (uint32_t i = 0; i < N_faces; i++) {
		(*faces)[i] = (void*) memblock;
		memblock += sizeof(face_t);
	}
}

static void init_global_matrix(vcn_sparse_t **K, const nb_graph_t *trg_x_vol,
			       const nb_partition_t *intmsh)
{
	uint32_t memsize = nb_graph_get_memsize();
	nb_graph_t *graph = NB_SOFT_MALLOC(memsize);

	nb_graph_init(graph);
	nb_cvfa_get_adj_graph(intmsh, trg_x_vol, graph);
	nb_graph_force_symmetry(graph);

	*K = vcn_sparse_create(graph, NULL, 2);

	nb_graph_finish(graph);
	NB_SOFT_FREE(memsize, graph);
}

static void load_faces(const nb_partition_t *part,
		       const nb_partition_t *intmsh,
		       const nb_graph_t *trg_x_vol,
		       face_t **faces)
{
	get_face_elems(part, faces);

	uint32_t N_faces = nb_partition_get_N_edges(part);
	for (uint32_t i = 0; i < N_faces; i++) {
		nb_partition_edge_get_normal(part, i, faces[i]->nf);
		
		uint32_t id1 = nb_partition_edge_get_1n(part, i);
		uint32_t id2 = nb_partition_edge_get_2n(part, i);
		faces[i]->x1[0] = nb_partition_node_get_x(part, id1);
		faces[i]->x1[1] = nb_partition_node_get_y(part, id1);
		faces[i]->x2[0] = nb_partition_node_get_x(part, id2);
		faces[i]->x2[1] = nb_partition_node_get_y(part, id2);

		load_subfaces(faces, i, intmsh, trg_x_vol);
	}
}

static void get_face_elems(const nb_partition_t *part, face_t **faces)
{
	uint32_t N_elems = nb_partition_get_N_elems(part);
	for (uint32_t i = 0; i < N_elems; i++) {
		uint16_t N_adj = nb_partition_elem_get_N_adj(part, i);
		for (uint16_t j = 0; j < N_adj; j++)
			define_face_elems(part, faces, i, j);
	}
}

static void define_face_elems(const nb_partition_t *part,
			      face_t **faces, uint32_t elem_id,
			      uint32_t local_face_id)
{
	uint32_t N_elems = nb_partition_get_N_elems(part);
	uint32_t ngb_id = nb_partition_elem_get_ngb(part, elem_id,
						    local_face_id);
	if (ngb_id >= N_elems) {
		uint32_t face_id = 
			nb_partition_elem_find_edge(part, elem_id,
						    local_face_id);
		faces[face_id]->elems[0] = elem_id;
		faces[face_id]->elems[1] = N_elems;
	} else if (elem_id < ngb_id) {
		uint32_t face_id = 
			nb_partition_elem_find_edge(part, elem_id,
						    local_face_id);
		uint32_t ev1 =
			nb_partition_elem_get_adj(part, elem_id,
						  local_face_id);
		uint32_t fv1 =
			nb_partition_edge_get_1n(part, face_id);
		if (ev1 == fv1) {
			faces[face_id]->elems[0] = elem_id;
			faces[face_id]->elems[1] = ngb_id;
		} else {
			faces[face_id]->elems[0] = ngb_id;
			faces[face_id]->elems[1] = elem_id;
		}
	}
}

static void load_subfaces(face_t **faces, uint32_t face_id,
			  const nb_partition_t *const intmsh,
			  const nb_graph_t *trg_x_vol)
{
	nb_container_type cnt_type = NB_QUEUE;
	uint32_t bank_size = nb_membank_get_memsize();
	uint32_t memsize = bank_size +
		nb_container_get_memsize(cnt_type);
	char *memblock = NB_SOFT_MALLOC(memsize);
	
	nb_membank_t *membank = (void*) memblock;
	nb_container_t *subfaces = (void*) (memblock + bank_size);

	nb_membank_init(membank, sizeof(subface_t));
	nb_container_init(subfaces, cnt_type);

	uint32_t elem_id = faces[face_id]->elems[0];
	uint16_t N_trg = trg_x_vol->N_adj[elem_id];
	uint32_t *trg_adj = trg_x_vol->adj[elem_id];

	uint8_t end_trg = 0;
	for (uint16_t i = 0; i < N_trg; i++) {
		uint8_t N_int = add_subface_if_intersected(membank,
							   intmsh, trg_adj,
							   faces, i, face_id,
							   subfaces);
		
		if (1 == N_int)
			end_trg += 1;
			
	}

	if (end_trg == 1)
		add_subface_in_closest_trg(membank, intmsh, faces,
					   face_id, subfaces);

	set_subfaces(membank, faces[face_id], subfaces);

	nb_container_finish(subfaces);
	nb_membank_finish(membank);
	NB_SOFT_FREE(memsize, memblock);
}

static uint8_t add_subface_if_intersected(nb_membank_t *membank,
					  const nb_partition_t *intmsh,
					  const uint32_t *trg_adj,
					  face_t **faces, uint32_t elem_trg_id,
					  uint32_t face_id,
					  nb_container_t *subfaces)
{
	face_t *face = faces[face_id];
	uint32_t trg_id = trg_adj[elem_trg_id];

	double t1[2], t2[2], t3[2];
	load_triangle_points(intmsh, trg_id, t1, t2, t3);

	uint8_t N_int = 0;
	double xp[4], p[2];

	if (vcn_utils2D_are_sgm_intersected(face->x1, face->x2, t1, t2, p)) {
		xp[N_int * 2] = p[0];
		xp[N_int*2+1] = p[1];
		N_int += 1;
	}

	if (vcn_utils2D_are_sgm_intersected(face->x1, face->x2, t2, t3, p)) {
		xp[N_int * 2] = p[0];
		xp[N_int*2+1] = p[1];
		N_int += 1;
	}
	
	if (2 > N_int) {
		if (vcn_utils2D_are_sgm_intersected(face->x1, face->x2,
						    t3, t1, p)) {
			xp[N_int * 2] = p[0];
			xp[N_int*2+1] = p[1];
			N_int += 1;
		}
	}

	if (1 == N_int) {
		if (vcn_utils2D_pnt_lies_in_trg(t1, t2, t3, face->x1)) {
			xp[2] = face->x1[0];
			xp[3] = face->x1[1];
		} else {
			xp[2] = face->x2[0];
			xp[3] = face->x2[1];			
		}
	}
	
	if (0 < N_int) {
		subface_t *subface = nb_membank_allocate_mem(membank);
		subface->N_int = N_int;
		memcpy(subface->x1, xp, 2 * sizeof(*xp));
		memcpy(subface->x2, &(xp[2]), 2 * sizeof(*xp));
		subface->trg_id = trg_id;
		nb_container_insert(subfaces, subface);
	}
	
	return N_int;
}

static void add_subface_in_closest_trg(nb_membank_t *membank,
				       const nb_partition_t *intmsh,
				       face_t **faces, uint32_t face_id,
				       nb_container_t *subfaces)
{
	face_t *face = faces[face_id];

	double alone[2];
	get_face_vtx_outside_intmsh(subfaces, face, alone);

	double p[2];
	uint16_t closest_id =
		get_face_closest_intersection_to_intmsh(subfaces,
							alone, p);

	subface_t *subface = nb_membank_allocate_mem(membank);
	subface->N_int = 1;
	memcpy(subface->x1, p, 2 * sizeof(*p));
	memcpy(subface->x2, alone, 2 * sizeof(*alone));
	subface->trg_id = closest_id;
	nb_container_insert(subfaces, subface);
}

static void get_face_vtx_outside_intmsh(const nb_container_t *subfaces,
					const face_t *face,
					double alone[2])
{
	bool s1_is_outside = true;

	nb_iterator_t * iter = alloca(nb_iterator_get_memsize());
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, subfaces);
	while (nb_iterator_has_more(iter)) {
		const subface_t *subface = nb_iterator_get_next(iter);
		if (1 == subface->N_int) {
			double e1 = fabs(face->x1[0] - subface->x2[0]);
			double e2 = fabs(face->x1[1] - subface->x2[1]);
			if (e1 < 1e-20 && e2 < 1e-20)
				s1_is_outside = false;
			break;
		}
	}
	nb_iterator_finish(iter);

	if (s1_is_outside)
		memcpy(alone, face->x1, 2 * sizeof(*alone));
	else
		memcpy(alone, face->x2, 2 * sizeof(*alone));
}

static uint16_t get_face_closest_intersection_to_intmsh
					(const nb_container_t *subfaces,
					 const double alone[2],
					 double p[2])
/* PENDING: SLOW FUNCTION (CALCULATE DIST AND DUPLICATE CHECKS) */
{
	uint16_t sf_id = 9999;
	double min = 1e30;
	nb_iterator_t * iter = alloca(nb_iterator_get_memsize());
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, subfaces);
	while (nb_iterator_has_more(iter)) {
		const subface_t *subface = nb_iterator_get_next(iter);
		double d = vcn_utils2D_get_dist(alone, subface->x1);
		if (d < min) {
			min = d;
			memcpy(p, subface->x1, 2 * sizeof(*p));
			sf_id = subface->trg_id;
		}
		d = vcn_utils2D_get_dist(alone, subface->x2);
		if (d < min) {
			min = d;
			memcpy(p, subface->x2, 2 * sizeof(*p));
			sf_id = subface->trg_id;
		}
	}
	nb_iterator_finish(iter);
	return sf_id;
}

static void set_subfaces(nb_membank_t *membank, face_t *face,
			 nb_container_t *subfaces)
{
	uint32_t N_sf = nb_container_get_length(subfaces);
	face->N_sf = N_sf;
	if (N_sf > 0) {
		uint32_t subfaces_size = N_sf * (sizeof(void*) +
						 sizeof(subface_t));
		char *memblock = nb_allocate_mem(subfaces_size);
		face->subfaces = (void*) memblock;
		memblock += N_sf * sizeof(void*);

		uint16_t cnt = 0;
		while (nb_container_is_not_empty(subfaces)) {
			subface_t *aux_sf =
				nb_container_delete_first(subfaces);
			subface_t *sf = (void*) memblock;
			memblock += sizeof(subface_t);

			sf->N_int = aux_sf->N_int;
			sf->trg_id = aux_sf->trg_id;
			memcpy(sf->x1, aux_sf->x1, 2 * sizeof(*(sf->x1)));
			memcpy(sf->x2, aux_sf->x2, 2 * sizeof(*(sf->x2)));

			face->subfaces[cnt] = sf;
			cnt += 1;

			nb_membank_free_mem(membank, aux_sf);
		}
	} else {
		face->subfaces = NULL;
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
				      const nb_partition_t *intmsh,
				      face_t **faces,
				      const nb_material_t *material,
				      nb_analysis2D_t analysis2D,
				      nb_analysis2D_params *params2D)
{
	vcn_sparse_reset(K);
	uint32_t N_faces = nb_partition_get_N_edges(part);
	for (uint32_t i = 0; i < N_faces; i++) {
		assemble_face(K, part, intmsh, faces[i], material,
			      analysis2D, params2D);
	}
}

static void assemble_face(vcn_sparse_t *K,
			  const nb_partition_t *const part,
			  const nb_partition_t *intmsh,
			  face_t *face,
			  const nb_material_t *material,
			  nb_analysis2D_t analysis2D,
			  nb_analysis2D_params *params2D)
{	
	double D[4];
	nb_pde_get_constitutive_matrix(D, material, analysis2D);

	uint32_t N_elems = nb_partition_get_N_elems(part);
	if (face->elems[1] < N_elems)
		assemble_internal_face(K, part, intmsh, face,
				       D, params2D);
}

static void assemble_internal_face(vcn_sparse_t *K,
				   const nb_partition_t *const part,
				   const nb_partition_t *intmsh,
				   face_t *face, const double D[4],
				   nb_analysis2D_params *params2D)
{
	if (0 && 0 < face->N_sf)/* TEMPORAL */
		integrate_subfaces(K, part, intmsh, face, D, params2D);
	else
		integrate_pairwise(K, part, face, D, params2D);
}

static void integrate_subfaces(vcn_sparse_t *K,
			       const nb_partition_t *const part,
			       const nb_partition_t *intmsh,
			       face_t *face, const double D[4],
			       nb_analysis2D_params *params2D)
{
	uint32_t memsize = 12 * sizeof(double); /* (N = 3) * 4 */
	char* memblock = NB_SOFT_MALLOC(memsize);
	double *Kf = (void*) memblock;
	
	uint16_t N_sf = face->N_sf;
	for (uint16_t i = 0; i < N_sf; i++) {
		integrate_Kf(part, intmsh, face, i, D, params2D, Kf);
		add_Kf_to_K(face, intmsh, i, Kf, K);
	}

	NB_SOFT_FREE(memsize, memblock);
}

static void integrate_Kf(const nb_partition_t *const part,
			 const nb_partition_t *intmsh, face_t *face,
			 uint16_t subface_id, const double D[4],
			 nb_analysis2D_params *params2D, double *Kf)
{
	memset(Kf, 0, 12 * sizeof(double));
	double t1[2], t2[2], t3[2];
	load_triangle_points(intmsh, face->subfaces[subface_id]->trg_id,
			     t1, t2, t3);

	double iJ[4];
	double detJ = subface_get_inverse_jacobian(t1, t2, t3, iJ);

	subface_t *subface = face->subfaces[subface_id];
	double lfn = subface_get_normalized_length(subface, t1, t2, t3);

	double *nf = face->nf;

	double factor = lfn * detJ * params2D->thickness;
	for (uint8_t i = 0; i < 3; i++) {
		double grad_xi[2];
		subface_get_normalized_grad(i, grad_xi);
		double grad[2];
		subface_get_grad(iJ, grad_xi, grad);
		double Kfi[4];
		subface_get_nodal_contribution(D, nf, grad, Kfi);
		Kf[i * 2] += factor * Kfi[0];
		Kf[i*2+1] += factor * Kfi[1];
		Kf[6 + i * 2] += factor * Kfi[2];
		Kf[6 + i*2+1] += factor * Kfi[3];
	}
}

static void load_triangle_points(const nb_partition_t *intmsh,
				 uint32_t trg_id, double t1[2],
				 double t2[2], double t3[2])
{
	uint32_t id1 = nb_partition_elem_get_adj(intmsh, trg_id, 0);
	uint32_t id2 = nb_partition_elem_get_adj(intmsh, trg_id, 1);
	uint32_t id3 = nb_partition_elem_get_adj(intmsh, trg_id, 2);

	t1[0] = nb_partition_node_get_x(intmsh, id1);
	t1[1] = nb_partition_node_get_y(intmsh, id1);

	t2[0] = nb_partition_node_get_x(intmsh, id2);
	t2[1] = nb_partition_node_get_y(intmsh, id2);

	t3[0] = nb_partition_node_get_x(intmsh, id3);
	t3[1] = nb_partition_node_get_y(intmsh, id3);
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
	get_normalized_point(t1, t2, t3, subface->x1, psi1);

	double psi2[2];
	get_normalized_point(t1, t2, t3, subface->x2, psi2);

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

static void add_Kf_to_K(face_t *face, const nb_partition_t *intmsh,
			uint16_t subface_id, const double *Kf,
			vcn_sparse_t *K)
{
	uint32_t i = face->elems[0];
	uint32_t j = face->elems[1];
	uint32_t trg_id = face->subfaces[subface_id]->trg_id;
	for (uint32_t m = 0; m < 3; m++) {
		uint32_t k = nb_partition_elem_get_adj(intmsh, trg_id, m);
		vcn_sparse_add(K, i * 2, k * 2, -Kf[m * 2]);
		vcn_sparse_add(K, i * 2, k*2+1, -Kf[m*2+1]);
		vcn_sparse_add(K, i*2+1, k * 2, -Kf[6 + m * 2]);
		vcn_sparse_add(K, i*2+1, k*2+1, -Kf[6 + m*2+1]);

		vcn_sparse_add(K, j * 2, k * 2, Kf[m * 2]);
		vcn_sparse_add(K, j * 2, k*2+1, Kf[m*2+1]);
		vcn_sparse_add(K, j*2+1, k * 2, Kf[6 + m * 2]);
		vcn_sparse_add(K, j*2+1, k*2+1, Kf[6 + m*2+1]);
	}
}

static void integrate_pairwise(vcn_sparse_t *K,
			       const nb_partition_t *part,
			       face_t *face, const double D[4],
			       nb_analysis2D_params *params2D)
{
	uint32_t memsize = 8 * sizeof(double); /* (N = 2) * 4 */
	char* memblock = NB_SOFT_MALLOC(memsize);
	double *Kf = (void*) memblock;

	integrate_Kf_pairwise(part, face, D, params2D, Kf);
	add_Kf_to_K_pairwise(face, Kf, K);

	NB_SOFT_FREE(memsize, memblock);
}

static void integrate_Kf_pairwise(const nb_partition_t *const part,
				  face_t *face, const double D[4],
				  nb_analysis2D_params *params2D, double *Kf)
{
	double c1[2], c2[2];
	c1[0] = nb_partition_elem_get_x(part, face->elems[0]);
	c1[1] = nb_partition_elem_get_y(part, face->elems[0]);
	c2[0] = nb_partition_elem_get_x(part, face->elems[1]);
	c2[1] = nb_partition_elem_get_y(part, face->elems[1]);

	double *nf = face->nf;

	double lf = vcn_utils2D_get_dist(face->x1, face->x2);
	double factor = lf * params2D->thickness;
	for (uint8_t i = 0; i < 2; i++) {
		double grad[2];
		if (0 == i)
			face_get_grad_pairwise(c1, c2, grad);
		else
			face_get_grad_pairwise(c2, c1, grad);
		double Kfi[4];
		subface_get_nodal_contribution(D, nf, grad, Kfi);
		Kf[i * 2] += factor * Kfi[0];
		Kf[i*2+1] += factor * Kfi[1];
		Kf[4 + i * 2] += factor * Kfi[2];
		Kf[4 + i*2+1] += factor * Kfi[3];
	}
}

static void face_get_grad_pairwise(const double c1[2], const double c2[2],
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

static void add_Kf_to_K_pairwise(face_t *face, const double *Kf,
				 vcn_sparse_t *K)
{
	uint32_t i = face->elems[0];
	uint32_t j = face->elems[1];
	for (uint32_t m = 0; m < 2; m++) {
		uint32_t k = face->elems[m];
		vcn_sparse_add(K, i * 2, k * 2, -Kf[m * 2]);
		vcn_sparse_add(K, i * 2, k*2+1, -Kf[m*2+1]);
		vcn_sparse_add(K, i*2+1, k * 2, -Kf[4 + m * 2]);
		vcn_sparse_add(K, i*2+1, k*2+1, -Kf[4 + m*2+1]);

		vcn_sparse_add(K, j * 2, k * 2, Kf[m * 2]);
		vcn_sparse_add(K, j * 2, k*2+1, Kf[m*2+1]);
		vcn_sparse_add(K, j*2+1, k * 2, Kf[4 + m * 2]);
		vcn_sparse_add(K, j*2+1, k*2+1, Kf[4 + m*2+1]);
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
			   face_t **faces,
			   const nb_partition_t *const part,
			   const nb_partition_t *intmsh,
			   const nb_bcond_t *const bcond,
			   const double *disp)
{
	uint32_t N_faces = nb_partition_get_N_edges(part);

 	for (uint32_t i = 0; i < N_faces; i++) {
		get_face_strain(faces, i, part, intmsh, bcond,
				disp, strain, boundary_mask);
	}
}

static void get_face_strain(face_t **faces, uint32_t face_id,
			    const nb_partition_t *const part,
			    const nb_partition_t *intmsh,
			    const nb_bcond_t *const bcond,
			    const double *disp,
			    double *strain,
			    char *boundary_mask)
{
	uint32_t N_elems = nb_partition_get_N_elems(part);
	if (faces[face_id]->elems[1] < N_elems) {
		boundary_mask[face_id] = 0;
		get_internal_face_strain(faces, face_id, part, intmsh,
					 disp, strain);
	} else {
		boundary_mask[face_id] = 1;
		get_boundary_face_strain(faces, face_id, part,
					 bcond, disp, strain);
	}
}

static void get_internal_face_strain(face_t **faces, uint32_t face_id,
				     const nb_partition_t *const part,
				     const nb_partition_t *intmsh,
				     const double *disp, double *strain)
{
	if (0 < faces[face_id]->N_sf)
		get_subfaces_strain(faces, face_id, part, intmsh,
				    disp, strain);
	else
		get_pairwise_strain(faces, face_id, part, disp, strain);
}

static void get_subfaces_strain(face_t **faces, uint32_t face_id,
				const nb_partition_t *const part,
				const nb_partition_t *intmsh,
				const double *disp, double *strain)
{
	
	memset(&(strain[face_id*3]), 0, 3 * sizeof(*strain));
	face_t *face = faces[face_id];
	for (uint16_t i = 0; i < face->N_sf; i++) {
		subface_t *subface = face->subfaces[i];
		subface_sum_strain_in_trg(part, intmsh, face_id,
					  subface, disp, strain);
	}
	double length = vcn_utils2D_get_dist(face->x1, face->x2);
	strain[face_id * 3] /= length;
	strain[face_id*3+1] /= length;
	strain[face_id*3+2] /= length;
}

static void subface_sum_strain_in_trg(const nb_partition_t *const part,
				      const nb_partition_t *intmsh,
				      uint32_t face_id,
				      const subface_t *subface,
				      const double *disp, double *strain)
{
	double t1[2], t2[2], t3[2];
	load_triangle_points(intmsh, subface->trg_id, t1, t2, t3);

	double iJ[4];
	double detJ = subface_get_inverse_jacobian(t1, t2, t3, iJ);

	double lfn = subface_get_normalized_length(subface, t1, t2, t3);

	double factor = lfn * detJ;
	for (uint8_t i = 0; i < 3; i++) {
		double grad_xi[2];
		subface_get_normalized_grad(i, grad_xi);
		double grad[2];
		subface_get_grad(iJ, grad_xi, grad);
		strain[face_id * 3] += factor * (grad[0] * disp[i * 2]);
		strain[face_id*3+1] += factor * (grad[1] * disp[i*2+1]);
		strain[face_id*3+2] += factor * (grad[1] * disp[i * 2] +
						 grad[0] * disp[i*2+1]);
	}
}

static void get_pairwise_strain(face_t **faces, uint32_t face_id,
				const nb_partition_t *part,
				const double *disp, double *strain)
{
	face_t *face = faces[face_id];

	double c1[2], c2[2];
	c1[0] = nb_partition_elem_get_x(part, face->elems[0]);
	c1[1] = nb_partition_elem_get_y(part, face->elems[0]);
	c2[0] = nb_partition_elem_get_x(part, face->elems[1]);
	c2[1] = nb_partition_elem_get_y(part, face->elems[1]);

	double *nf = face->nf;

	double lf = vcn_utils2D_get_dist(face->x1, face->x2);

	for (uint8_t i = 0; i < 2; i++) {
		double grad[2];
		if (0 == i)
			face_get_grad_pairwise(c1, c2, grad);
		else
			face_get_grad_pairwise(c2, c1, grad);
		uint32_t id = face->elems[i];
		strain[face_id * 3] += lf * (grad[0] * disp[id * 2]);
		strain[face_id*3+1] += lf * (grad[1] * disp[id*2+1]);
		strain[face_id*3+2] += lf * (grad[1] * disp[id * 2] +
					     grad[0] * disp[id*2+1]);
	}
}

static void get_boundary_face_strain(face_t **faces, uint32_t face_id,
				     const nb_partition_t *const part,
				     const nb_bcond_t *bcond,
				     const double *disp, double *strain)
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

static void finish_faces(uint32_t N_faces, face_t **faces)
{
	for (uint32_t i = 0; i < N_faces; i++) {
		nb_free_mem(faces[i]->subfaces);
	}
}
