#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "nb/math_bot.h"
#include "nb/memory_bot.h"
#include "nb/container_bot.h"
#include "nb/solver_bot.h"
#include "nb/geometric_bot.h"
#include "nb/graph_bot.h"
#include "nb/pde_bot.h"

#include "../calculation_points.h"
#include "../integration_mesh.h"

#include "elasticity2D.h"
#include "set_bconditions.h"

#define SMOOTH 2

static uint32_t get_cvfa_memsize(uint32_t N_elems, uint32_t N_faces);
static void distribute_cvfa_memory(char *memblock, uint32_t N_elems,
				   uint32_t N_faces, double **xc, double **F,
				   nb_mesh2D_t **intmsh, nb_graph_t **trg_x_vol,
				   face_t ***faces, nb_glquadrature_t *glq);
static void init_global_matrix(nb_sparse_t **K, const nb_graph_t *trg_x_vol,
			       const nb_mesh2D_t *intmsh);
static void load_faces(const nb_mesh2D_t *mesh,
		       const nb_mesh2D_t *intmsh,
		       const nb_graph_t *trg_x_vol,
		       face_t **faces);
static void get_face_elems(const nb_mesh2D_t *mesh, face_t **faces);
static void define_face_elems(const nb_mesh2D_t *mesh,
			      face_t **faces, uint32_t elem_id,
			      uint16_t local_face_id);
static void load_subfaces(face_t **faces, uint32_t face_id,
			  const nb_mesh2D_t *const intmsh,
			  const nb_graph_t *trg_x_vol);
static uint8_t add_subface_if_intersected(nb_membank_t *membank,
					  const nb_mesh2D_t *intmsh,
					  const uint32_t *trg_adj,
					  face_t **faces, uint32_t elem_trg_id,
					  uint32_t face_id,
					  nb_container_t *subfaces);

static bool is_subface_inside_trg(const nb_mesh2D_t *intmsh,
				  uint16_t N_trg, const uint32_t *trg_adj,
				  face_t **faces, uint32_t face_id,
				  uint32_t *trg_id);
static void add_subface_inside_trg(nb_membank_t *membank, face_t **faces,
				   uint32_t face_id, uint32_t trg_id,
				   nb_container_t *subfaces);
static void add_subface_outside_trg(nb_membank_t *membank, face_t **faces,
				    uint32_t face_id, nb_container_t *subfaces);
static void add_subfaces_pairwise_ends(nb_membank_t *membank,
				       face_t **faces, uint32_t face_id,
				       nb_container_t *subfaces);
static void add_subface_pairwise(nb_membank_t *membank,
				       face_t **faces, uint32_t face_id,
				       nb_container_t *subfaces);

static void get_face_vtx_outside_intmsh(const nb_container_t *subfaces,
					const face_t *face,
					double alone[2]);
static uint32_t get_face_closest_intersection_to_intmsh
					(const nb_container_t *subfaces,
					 const double alone[2],
					 double p[2]);
static void set_subfaces(nb_membank_t *membank, face_t *face,
			 nb_container_t *subfaces);
static void finish_faces(uint32_t N_faces, face_t **faces);

int nb_cvfa_compute_2D_Solid_Mechanics
			(const nb_mesh2D_t *const mesh,
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
	uint32_t N_elems = nb_mesh2D_get_N_elems(mesh);
	uint32_t N_faces = nb_mesh2D_get_N_edges(mesh);
	uint32_t memsize = get_cvfa_memsize(N_elems, N_faces);
	char *memblock = nb_soft_allocate_mem(memsize);
	double *xc;
	double *F;
	nb_mesh2D_t *intmsh;
	nb_graph_t *trg_x_vol;
	face_t **faces;
	nb_glquadrature_t glq;
	distribute_cvfa_memory(memblock, N_elems, N_faces, &xc, &F,
			       &intmsh, &trg_x_vol, &faces, &glq);

	nb_glquadrature_load(&glq, SMOOTH + 1);

  	nb_cvfa_set_calculation_points(mesh, xc);
	nb_cvfa_init_integration_mesh(intmsh);
	nb_cvfa_load_integration_mesh(intmsh, N_elems, xc);

	nb_graph_init(trg_x_vol);
	nb_cvfa_correlate_mesh_and_integration_mesh(mesh, intmsh,
						    trg_x_vol);

	nb_sparse_t *K;
	init_global_matrix(&K, trg_x_vol, intmsh);

	load_faces(mesh, intmsh, trg_x_vol, faces);

	status = nb_cvfa_solve_elasticity_equation(mesh, material, bcond,
						   SMOOTH, enable_self_weight,
						   gravity, analysis2D,
						   params2D, displacement,
						   intmsh, xc, faces, F,
						   K, &glq);
	if (status != 0)
		goto CLEAN_AND_EXIT;

	nb_cvfa_compute_strain(strain, boundary_mask, faces, mesh, SMOOTH,
			       intmsh, xc, bcond, displacement, &glq);

	status = 0;
CLEAN_AND_EXIT:
	finish_faces(N_faces, faces);
	nb_sparse_destroy(K);
	nb_graph_finish(trg_x_vol);
	nb_mesh2D_finish(intmsh);
	nb_soft_free_mem(memsize, memblock);
	return status;
}

static uint32_t get_cvfa_memsize(uint32_t N_elems, uint32_t N_faces)
{
	uint32_t system_size = 4 * N_elems * sizeof(double);
	uint32_t intmsh_size = nb_cvfa_get_integration_mesh_memsize();
	uint32_t graph_size = nb_graph_get_memsize();
	uint16_t Nq = SMOOTH + 1;
	uint32_t glq_size = 2 * Nq * sizeof(double);
	uint32_t faces_size = N_faces * (sizeof(void*) + sizeof(face_t));
	return graph_size + system_size + intmsh_size + faces_size + glq_size;
}

static void distribute_cvfa_memory(char *memblock, uint32_t N_elems,
				   uint32_t N_faces, double **xc, double **F,
				   nb_mesh2D_t **intmsh, nb_graph_t **trg_x_vol,
				   face_t ***faces, nb_glquadrature_t *glq)
{
	uint32_t system_size = 2 * N_elems * sizeof(double);
	uint32_t intmsh_size = nb_cvfa_get_integration_mesh_memsize();
	uint32_t graph_size = nb_graph_get_memsize();
	uint16_t Nq = SMOOTH + 1;
	uint32_t glq_size = 2 * Nq * sizeof(double);
	*F = (void*) memblock;
	*xc = (void*) (memblock + system_size);
	*intmsh = (void*) (memblock + 2 * system_size);
	*trg_x_vol = (void*) (memblock + 2 * system_size + intmsh_size);
	glq->x = (void*) (memblock + 2 * system_size +
			  intmsh_size + graph_size);
	glq->w = (void*) (memblock + 2 * system_size +
			  intmsh_size + graph_size + Nq * sizeof(double));
	*faces = (void*) (memblock + 2 * system_size +
			  intmsh_size + graph_size + glq_size);
	memblock +=  2 * system_size + intmsh_size + graph_size +
		glq_size + N_faces * sizeof(void*);
	for (uint32_t i = 0; i < N_faces; i++) {
		(*faces)[i] = (void*) (memblock + i * sizeof(face_t));
		memset((*faces)[i], 0, sizeof(face_t));
	}
}

static void init_global_matrix(nb_sparse_t **K, const nb_graph_t *trg_x_vol,
			       const nb_mesh2D_t *intmsh)
{
	uint32_t memsize = nb_graph_get_memsize();
	nb_graph_t *graph = nb_soft_allocate_mem(memsize);

	nb_graph_init(graph);
	nb_cvfa_get_adj_graph(intmsh, trg_x_vol, graph);

	/* Forces symmetry for LU decomposition */
	nb_graph_force_symmetry(graph);

	*K = nb_sparse_create(graph, NULL, 2);

	nb_graph_finish(graph);
	nb_soft_free_mem(memsize, graph);
}

static void load_faces(const nb_mesh2D_t *mesh,
		       const nb_mesh2D_t *intmsh,
		       const nb_graph_t *trg_x_vol,
		       face_t **faces)
{
	get_face_elems(mesh, faces);

	uint32_t N_faces = nb_mesh2D_get_N_edges(mesh);
	for (uint32_t i = 0; i < N_faces; i++) {
		nb_mesh2D_edge_get_normal(mesh, i, faces[i]->nf);		
		uint32_t id1 = nb_mesh2D_edge_get_1n(mesh, i);
		uint32_t id2 = nb_mesh2D_edge_get_2n(mesh, i);
		faces[i]->x1[0] = nb_mesh2D_node_get_x(mesh, id1);
		faces[i]->x1[1] = nb_mesh2D_node_get_y(mesh, id1);
		faces[i]->x2[0] = nb_mesh2D_node_get_x(mesh, id2);
		faces[i]->x2[1] = nb_mesh2D_node_get_y(mesh, id2);

		load_subfaces(faces, i, intmsh, trg_x_vol);
	}
}

static void get_face_elems(const nb_mesh2D_t *mesh, face_t **faces)
{
	uint32_t N_elems = nb_mesh2D_get_N_elems(mesh);
	for (uint32_t i = 0; i < N_elems; i++) {
		uint16_t N_adj = nb_mesh2D_elem_get_N_adj(mesh, i);
		for (uint16_t j = 0; j < N_adj; j++)
			define_face_elems(mesh, faces, i, j);
	}
}

static void define_face_elems(const nb_mesh2D_t *mesh,
			      face_t **faces, uint32_t elem_id,
			      uint16_t local_face_id)
{
	uint32_t N_elems = nb_mesh2D_get_N_elems(mesh);
	uint32_t ngb_id = nb_mesh2D_elem_get_ngb(mesh, elem_id,
						 local_face_id);
	if (ngb_id >= N_elems) {
		uint32_t face_id = 
			nb_mesh2D_elem_find_edge(mesh, elem_id,
						 local_face_id);
		faces[face_id]->elems[0] = elem_id;
		faces[face_id]->elems[1] = N_elems;
	} else if (elem_id < ngb_id) {
		uint32_t face_id = 
			nb_mesh2D_elem_find_edge(mesh, elem_id,
						 local_face_id);
		uint32_t ev1 =
			nb_mesh2D_elem_get_adj(mesh, elem_id,
					       local_face_id);
		uint32_t fv1 =
			nb_mesh2D_edge_get_1n(mesh, face_id);
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
			  const nb_mesh2D_t *const intmsh,
			  const nb_graph_t *trg_x_vol)
{
	nb_container_type cnt_type = NB_QUEUE;
	uint32_t bank_size = nb_membank_get_memsize();
	uint32_t memsize = bank_size +
		nb_container_get_memsize(cnt_type);
	char *memblock = nb_soft_allocate_mem(memsize);
	
	nb_membank_t *membank = (void*) memblock;
	nb_container_t *subfaces = (void*) (memblock + bank_size);

	nb_membank_init(membank, sizeof(subface_t));
	nb_container_init(subfaces, cnt_type);

	uint32_t elem_id = faces[face_id]->elems[0];
	uint16_t N_trg = trg_x_vol->N_adj[elem_id];
	uint32_t *trg_adj = trg_x_vol->adj[elem_id];

	uint8_t end_trg = 0;
	for (uint16_t i = 0; i < N_trg; i++) {
		uint8_t N_int = add_subface_if_intersected(membank, intmsh,
							   trg_adj, faces,
							   i, face_id,
							   subfaces);
		
		if (1 == N_int)
			end_trg += 1;
			
	}

	if (0 == end_trg) {
		uint8_t N_sf = nb_container_get_length(subfaces);
		if (0 == N_sf) {
			uint32_t trg_id;
			bool inside = is_subface_inside_trg(intmsh, N_trg,
							    trg_adj, faces,
							    face_id, &trg_id);
			if (inside)
				add_subface_inside_trg(membank, faces, face_id,
						       trg_id, subfaces);

			else
				add_subface_outside_trg(membank, faces, face_id,
							subfaces);
		} else {
			add_subfaces_pairwise_ends(membank, faces, face_id,
						   subfaces);
		}
	}

	if (1 == end_trg)
		add_subface_pairwise(membank, faces, face_id, subfaces);

	set_subfaces(membank, faces[face_id], subfaces);

	nb_container_finish(subfaces);
	nb_membank_finish(membank);
	nb_soft_free_mem(memsize, memblock);
}

static uint8_t add_subface_if_intersected(nb_membank_t *membank,
					  const nb_mesh2D_t *intmsh,
					  const uint32_t *trg_adj,
					  face_t **faces, uint32_t elem_trg_id,
					  uint32_t face_id,
					  nb_container_t *subfaces)
{
	/*   o---------o---------o
	 *    \       / \       /
	 *     \  +--/---\--+  /   Find face intersections.
	 *      \   /     \   /
	 *       \ /       \ /
	 *        o --------o
	 */
	face_t *face = faces[face_id];
	uint32_t trg_id = trg_adj[elem_trg_id];

	double t1[2], t2[2], t3[2];
	nb_cvfa_load_trg_points(intmsh, trg_id, t1, t2, t3);

	uint8_t N_int = 0;
	double xp[4], p[2];

	if (nb_utils2D_are_sgm_intersected(face->x1, face->x2, t1, t2, p)) {
		xp[N_int * 2] = p[0];
		xp[N_int*2+1] = p[1];
		N_int += 1;
	}

	if (nb_utils2D_are_sgm_intersected(face->x1, face->x2, t2, t3, p)) {
		xp[N_int * 2] = p[0];
		xp[N_int*2+1] = p[1];
		N_int += 1;
	}
	
	if (2 > N_int) {
		if (nb_utils2D_are_sgm_intersected(face->x1, face->x2,
						    t3, t1, p)) {
			xp[N_int * 2] = p[0];
			xp[N_int*2+1] = p[1];
			N_int += 1;
		}
	}

	if (1 == N_int) {
		if (nb_utils2D_pnt_lies_in_trg(t1, t2, t3, face->x1)) {
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

static bool is_subface_inside_trg(const nb_mesh2D_t *intmsh,
				  uint16_t N_trg, const uint32_t *trg_adj,
				  face_t **faces, uint32_t face_id,
				  uint32_t *trg_id)
{
	face_t *face = faces[face_id];

	*trg_id = nb_mesh2D_get_N_elems(intmsh);
	bool inside = false;
	for (uint16_t i = 0; i < N_trg; i++) {
		uint32_t id = trg_adj[i];

		double t1[2], t2[2], t3[2];
		nb_cvfa_load_trg_points(intmsh, id, t1, t2, t3);

		bool s1_in = nb_utils2D_pnt_lies_in_trg(t1, t2, t3, face->x1);
		bool s2_in = nb_utils2D_pnt_lies_in_trg(t1, t2, t3, face->x2);

		if (s1_in && s2_in) {
			inside = true;
			*trg_id = id;
			break;
		}
	}
	return inside;
}

static void add_subface_inside_trg(nb_membank_t *membank, face_t **faces,
				   uint32_t face_id, uint32_t trg_id,
				   nb_container_t *subfaces)
{
	/*   o---------o
	 *    \ +---+ /
	 *     \     /   Face contained in trg.
	 *      \   /
	 *       \ /
	 *        o
	 */
	face_t *face = faces[face_id];

	subface_t *subface = nb_membank_allocate_mem(membank);
	subface->N_int = 1;
	memcpy(subface->x1, face->x1, 2 * sizeof(*(face->x1)) );
	memcpy(subface->x2, face->x2, 2 * sizeof(*(face->x2)) );
	subface->trg_id = trg_id;
	nb_container_insert(subfaces, subface);
}

static void add_subface_outside_trg(nb_membank_t *membank, face_t **faces,
				   uint32_t face_id, nb_container_t *subfaces)
{
	/*   o---------o
	 *    \       /
	 *     \     /   Face not contained in
	 *      \   /    any trg
	 * +---+ \ / 
	 *        o
	 */
	face_t *face = faces[face_id];

	subface_t *subface = nb_membank_allocate_mem(membank);
	subface->N_int = 0;
	memcpy(subface->x1, face->x1, 2 * sizeof(*(face->x1)) );
	memcpy(subface->x2, face->x2, 2 * sizeof(*(face->x2)) );
	subface->trg_id = 0;
	nb_container_insert(subfaces, subface);
}

static void add_subfaces_pairwise_ends(nb_membank_t *membank,
				       face_t **faces, uint32_t face_id,
				       nb_container_t *subfaces)
{
	/*   o---------o
	 *    \       /
	 *     \     /     Face intersected by trg with
	 *    +-\---/-+    both ending points outside.
	 *       \ /
	 *        o
	 */
	add_subface_pairwise(membank, faces, face_id, subfaces);
	add_subface_pairwise(membank, faces, face_id, subfaces);

}

static void add_subface_pairwise(nb_membank_t *membank,
				 face_t **faces, uint32_t face_id,
				 nb_container_t *subfaces)
{
	/*   o---------o
	 *    \       /
	 *     \     /     Face intersected by trg with
	 *    +-\-+ /      an ending point outside.
	 *       \ /
	 *        o
	 */
	face_t *face = faces[face_id];

	double alone[2];
	get_face_vtx_outside_intmsh(subfaces, face, alone);

	double p[2];
	uint16_t closest_id =
		get_face_closest_intersection_to_intmsh(subfaces,
							alone, p);

	subface_t *subface = nb_membank_allocate_mem(membank);
	subface->N_int = 0;
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

	nb_iterator_t * iter = nb_allocate_on_stack(nb_iterator_get_memsize());
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

static uint32_t get_face_closest_intersection_to_intmsh
					(const nb_container_t *subfaces,
					 const double alone[2],
					 double p[2])
/* PENDING: SLOW FUNCTION (CALCULATE DIST AND DUPLICATE CHECKS) */
{
	uint32_t sf_id = 9999;
	double min = 1e30;
	nb_iterator_t * iter = nb_allocate_on_stack(nb_iterator_get_memsize());
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, subfaces);
	while (nb_iterator_has_more(iter)) {
		const subface_t *subface = nb_iterator_get_next(iter);
		double d = nb_utils2D_get_dist(alone, subface->x1);
		if (d < min) {
			min = d;
			memcpy(p, subface->x1, 2 * sizeof(*p));
			sf_id = subface->trg_id;
		}
		d = nb_utils2D_get_dist(alone, subface->x2);
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

void nb_cvfa_compute_stress_from_strain(const nb_mesh2D_t *mesh,
					const nb_material_t *const material,
					nb_analysis2D_t analysis2D,
					const double* strain,
					double* stress /* Output */)
{
	uint32_t N_faces = nb_mesh2D_get_N_edges(mesh);
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
