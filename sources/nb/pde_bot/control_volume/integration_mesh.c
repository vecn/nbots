#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "nb/memory_bot.h"
#include "nb/container_bot.h"
#include "nb/solver_bot.h"
#include "nb/geometric_bot.h"

#include "integration_mesh.h"

#define INTEGRATOR_TYPE NB_TRIAN

static void init_containers_trg_x_vol(nb_container_t **all_trg_x_vol,
				      nb_container_type cnt_type,
				      nb_membank_t *membank,
				      const nb_mesh2D_t *mesh,
				      const nb_mesh2D_t *intmsh);
static void vol_get_adj(const nb_mesh2D_t *mesh,
			const nb_mesh2D_t *intmsh,
			const nb_graph_t *trg_x_vtx,
			nb_membank_t *membank,
			uint32_t vol_id,
			nb_container_t *trg_adj);
static int8_t compare_ids(const void *ptr1, const void *ptr2);
static bool vol_intersects_trg(const nb_mesh2D_t *mesh,
			       const nb_mesh2D_t *intmsh,
			       uint32_t vol_id, uint32_t trg_id);
static void mesh_load_sgm_from_adj(const nb_mesh2D_t *mesh,
				   uint32_t elem_id, uint16_t adj_id,
				   double s1[2], double s2[2]);
static void put_neighbours_in_active(const nb_mesh2D_t *intmsh,
				     nb_membank_t *membank,
				     const nb_container_t *trg_adj,
				     const nb_container_t *out,
				     nb_container_t *active,
				     uint32_t id);
static void insert_in_active_if_dont_exists(nb_container_t *active,
					    nb_membank_t *membank,
					    const nb_container_t *trg_adj,
					    const nb_container_t *out,
					    uint32_t ngb);
static void trg_x_vol_allocate_adj(nb_graph_t *trg_x_vol,
				   nb_container_t **all_trg_x_vol);
static uint32_t trg_x_vol_get_N_adj(uint32_t N,
				    nb_container_t **all_trg_x_vol);
static void trg_x_vol_set_adj(nb_graph_t *trg_x_vol,
			      nb_container_t **all_trg_x_vol,
			      nb_membank_t *membank);
static void finish_containers_trg_x_vol(uint32_t N,
					nb_container_t **all_trg_x_vol,
					nb_membank_t *membank);
static void adj_graph_allocate_adj(nb_graph_t *graph,
				   const nb_graph_t *trg_x_vol,
				   const nb_mesh2D_t *intmsh);
static uint32_t adj_graph_get_N_adj(const nb_graph_t *trg_x_vol,
				    const nb_mesh2D_t *intmsh);
static void adj_graph_get_list_x_vol(const nb_graph_t *trg_x_vol,
				     const nb_mesh2D_t *intmsh,
				     uint32_t vol_id,
				     nb_membank_t *membank,
				     nb_container_t *list);
static void adj_graph_set_adj(nb_graph_t *graph,
			      const nb_graph_t *trg_x_vol,
			      const nb_mesh2D_t *intmsh);
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

uint32_t nb_cvfa_get_integration_mesh_memsize(void)
{
	return nb_mesh2D_get_memsize(INTEGRATOR_TYPE);
}

void nb_cvfa_init_integration_mesh(nb_mesh2D_t *intmsh)
{
	nb_mesh2D_init(intmsh, INTEGRATOR_TYPE);
}

void nb_cvfa_load_integration_mesh(nb_mesh2D_t *intmsh, uint32_t N,
				   const double *xc)
{
	uint32_t mesh_size = nb_tessellator2D_get_memsize();
	uint32_t perm_size = N * sizeof(uint32_t);
	uint32_t memsize = mesh_size + perm_size;
	char *memblock = nb_soft_allocate_mem(memsize);

	nb_tessellator2D_t *t2d = (void*) memblock;
	uint32_t *perm = (void*) (memblock + mesh_size);

	nb_tessellator2D_init(t2d);
	nb_tessellator2D_get_smallest_ns_alpha_complex(t2d, N, xc, 0.666);
	nb_mesh2D_load_from_tessellator2D(intmsh, t2d);
	nb_tessellator2D_finish(t2d);

	for (uint32_t i = 0; i < N; i++) {
		uint32_t id = nb_mesh2D_get_invtx(intmsh, i);
		perm[id] = i;
	}

	nb_mesh2D_set_nodal_permutation(intmsh, perm);

	nb_soft_free_mem(memsize, memblock);
}

void nb_cvfa_correlate_mesh_and_integration_mesh
					(const nb_mesh2D_t *mesh,
					 const nb_mesh2D_t *intmsh,
					 nb_graph_t *trg_x_vol)
{
	nb_container_type cnt_type = NB_SORTED;
	uint32_t N_elems = nb_mesh2D_get_N_elems(mesh);
	uint32_t cnt_size = nb_container_get_memsize(cnt_type);
	uint32_t bank_size = nb_membank_get_memsize();
	uint32_t memsize = N_elems * (cnt_size + sizeof(void*)) + bank_size;
	char *memblock = nb_soft_allocate_mem(memsize);
	nb_container_t **all_trg_x_vol = (void*) memblock;
	nb_membank_t *membank = (void*) (memblock + N_elems * 
					 (cnt_size + sizeof(void*)));

	nb_membank_init(membank, sizeof(uint32_t));

	init_containers_trg_x_vol(all_trg_x_vol, cnt_type, membank,
				  mesh, intmsh);
	
	trg_x_vol->N = N_elems;
	trg_x_vol_allocate_adj(trg_x_vol, all_trg_x_vol);
	trg_x_vol_set_adj(trg_x_vol, all_trg_x_vol, membank);

	finish_containers_trg_x_vol(N_elems, all_trg_x_vol, membank);

	nb_membank_finish(membank);
	nb_soft_free_mem(memsize, memblock);
}

static void init_containers_trg_x_vol(nb_container_t **all_trg_x_vol,
				      nb_container_type cnt_type,
				      nb_membank_t *membank,
				      const nb_mesh2D_t *mesh,
				      const nb_mesh2D_t *intmsh)
{
	uint32_t memsize = nb_graph_get_memsize();
	char *memblock = nb_soft_allocate_mem(memsize);
	nb_graph_t *trg_x_vtx = (void*) memblock;

	nb_graph_init(trg_x_vtx);
	nb_mesh2D_load_graph(intmsh, trg_x_vtx, NB_ELEMS_CONNECTED_TO_NODES);

	uint32_t N_elems = nb_mesh2D_get_N_elems(mesh);
	char *block = ((char*) all_trg_x_vol) + N_elems * sizeof(void*);
	uint32_t cnt_size = nb_container_get_memsize(cnt_type);
	for (uint32_t i = 0; i < N_elems; i++) {
		nb_container_t *trg_adj = (void*) block;

		all_trg_x_vol[i] = trg_adj;
		nb_container_init(trg_adj, cnt_type);

		vol_get_adj(mesh, intmsh, trg_x_vtx, membank, i, trg_adj);

		block += cnt_size;
	}

	nb_graph_finish(trg_x_vtx);

	nb_soft_free_mem(memsize, memblock);
}

static void vol_get_adj(const nb_mesh2D_t *mesh,
			const nb_mesh2D_t *intmsh,
			const nb_graph_t *trg_x_vtx,
			nb_membank_t *membank,
			uint32_t vol_id,
			nb_container_t *trg_adj)
{
	nb_container_type cnt_type = NB_SORTED;
	uint32_t cnt_size = nb_container_get_memsize(cnt_type);
	uint32_t active_size = nb_container_get_memsize(cnt_type);
	uint32_t memsize = cnt_size + active_size;
	char *memblock = nb_soft_allocate_mem(memsize);
	nb_container_t *out = (void*) memblock;
	nb_container_t *active = (void*) (memblock + cnt_size);

	nb_container_init(out, cnt_type);
	nb_container_init(active, cnt_type);

	nb_container_set_comparer(trg_adj, compare_ids);
	nb_container_set_comparer(out, compare_ids);
	nb_container_set_comparer(active, compare_ids);

	uint32_t *id = nb_membank_allocate_mem(membank);
	*id = trg_x_vtx->adj[vol_id][0];
	nb_container_insert(active, id);
	
	while (nb_container_is_not_empty(active)) {
		id = nb_container_delete_first(active);
		bool intersection = vol_intersects_trg(mesh, intmsh,
						       vol_id, *id);
		if (intersection) {
			nb_container_insert(trg_adj, id);
			put_neighbours_in_active(intmsh, membank, trg_adj, out,
						 active, *id);
		} else {
			nb_container_insert(out, id);
		}
	}

	while (nb_container_is_not_empty(out)) {
		id = nb_container_delete_first(out);
		nb_membank_free_mem(membank, id);
	}

	nb_container_finish(out);
	nb_container_finish(active);

	nb_soft_free_mem(memsize, memblock);
}

static int8_t compare_ids(const void *ptr1, const void *ptr2)
{
	const uint32_t *id1 = ptr1;
	const uint32_t *id2 = ptr2;
	int8_t out;
	if (*id1 < *id2)
		out = -1;
	else if (*id1 > *id2)
		out = 1;
	else
		out = 0;
	return out;
}

static bool vol_intersects_trg(const nb_mesh2D_t *mesh,
			       const nb_mesh2D_t *intmsh,
			       uint32_t vol_id, uint32_t trg_id)
{
	bool out = false;

	double a1[2], a2[2], b1[2], b2[2];
	
	uint16_t N_adj1 = nb_mesh2D_elem_get_N_adj(mesh, vol_id);
	for (uint16_t i = 0; i < N_adj1; i++) {
		mesh_load_sgm_from_adj(mesh, vol_id, i, a1, a2);
		uint16_t N_adj2 = nb_mesh2D_elem_get_N_adj(intmsh, trg_id);
		for (uint16_t j = 0; j < N_adj2; j++) {
			mesh_load_sgm_from_adj(intmsh, trg_id, j, b1, b2);
			
			nb_intersect_t status =
				nb_utils2D_get_sgm_intersection(a1, a2, b1,
								b2, NULL);
			if (NB_INTERSECTED == status ||
			    NB_PARALLEL == status) {
				out = true;
				goto EXIT;
			}
		}
	}
EXIT:
	return out;
}

static void mesh_load_sgm_from_adj(const nb_mesh2D_t *mesh,
				   uint32_t elem_id, uint16_t adj_id,
				   double s1[2], double s2[2])
{
	uint32_t N_adj = nb_mesh2D_elem_get_N_adj(mesh, elem_id);
	uint32_t id1 = nb_mesh2D_elem_get_adj(mesh, elem_id, adj_id);
	uint32_t id2 = nb_mesh2D_elem_get_adj(mesh, elem_id,
						 (adj_id + 1) % N_adj);
	s1[0] = nb_mesh2D_node_get_x(mesh, id1);
	s1[1] = nb_mesh2D_node_get_y(mesh, id1);

	s2[0] = nb_mesh2D_node_get_x(mesh, id2);
	s2[1] = nb_mesh2D_node_get_y(mesh, id2);
}

static void put_neighbours_in_active(const nb_mesh2D_t *intmsh,
				     nb_membank_t *membank,
				     const nb_container_t *trg_adj,
				     const nb_container_t *out,
				     nb_container_t *active, uint32_t id)
{
	uint16_t N_adj = nb_mesh2D_elem_get_N_adj(intmsh, id);
	for (uint32_t i = 0; i < N_adj; i++) {
		if (nb_mesh2D_elem_has_ngb(intmsh, id, i)) {
			uint32_t ngb = nb_mesh2D_elem_get_ngb(intmsh,
								 id, i);
			insert_in_active_if_dont_exists(active, membank,
							trg_adj, out, ngb);
		}
	}
}

static void insert_in_active_if_dont_exists(nb_container_t *active,
					    nb_membank_t *membank,
					    const nb_container_t *trg_adj,
					    const nb_container_t *out,
					    uint32_t ngb)
{
	if (NULL == nb_container_exist(trg_adj, &ngb)) {
		if (NULL == nb_container_exist(out, &ngb)) {
			if (NULL == nb_container_exist(active, &ngb)) {
				uint32_t *id = 
					nb_membank_allocate_mem(membank);
				*id = ngb;
				nb_container_insert(active, id);
			}
		}
	}
}

static void trg_x_vol_allocate_adj(nb_graph_t *trg_x_vol,
				   nb_container_t **all_trg_x_vol)
{
	uint32_t memsize_N_adj = trg_x_vol->N * sizeof(*(trg_x_vol->N_adj));
	uint32_t N_adj = trg_x_vol_get_N_adj(trg_x_vol->N, all_trg_x_vol);
	uint32_t memsize_adj = trg_x_vol->N * sizeof(*(trg_x_vol->adj)) +
		N_adj * sizeof(**(trg_x_vol->adj));
	char *memblock = nb_allocate_mem(memsize_N_adj + memsize_adj);
	trg_x_vol->N_adj = (void*) memblock;
	trg_x_vol->adj = (void*) (memblock + memsize_N_adj);
}

static uint32_t trg_x_vol_get_N_adj(uint32_t N,
				    nb_container_t **all_trg_x_vol)
{
	uint32_t N_adj = 0;
	for (uint32_t i = 0; i < N; i++) {
		nb_container_t *trg_adj = all_trg_x_vol[i];
		N_adj += nb_container_get_length(trg_adj);
	}
	return N_adj;	
}

static void trg_x_vol_set_adj(nb_graph_t *trg_x_vol,
			      nb_container_t **all_trg_x_vol,
			      nb_membank_t *membank)
{
	uint32_t mem_used = trg_x_vol->N * sizeof(*(trg_x_vol->N_adj)) +
		trg_x_vol->N * sizeof(*(trg_x_vol->adj));
	char *block = (char*) trg_x_vol->N_adj + mem_used;

	for (uint32_t i = 0; i < trg_x_vol->N; i++) {
		nb_container_t *trg_adj = all_trg_x_vol[i];
		uint16_t N_adj = nb_container_get_length(trg_adj);
		trg_x_vol->N_adj[i] = N_adj;

		trg_x_vol->adj[i] = (void*) block;
		block += N_adj * sizeof(**(trg_x_vol->adj));

		uint16_t j = 0;
		while (nb_container_is_not_empty(trg_adj)) {
			uint32_t *id = nb_container_delete_first(trg_adj);

			trg_x_vol->adj[i][j] = *id;
			j += 1;

			nb_membank_free_mem(membank, id);
		}
	}
}

static void finish_containers_trg_x_vol(uint32_t N,
					nb_container_t **all_trg_x_vol,
					nb_membank_t *membank)
{
	for (uint32_t i = 0; i < N; i++) {
		nb_container_t *trg_adj = all_trg_x_vol[i];
		nb_container_finish(trg_adj);
	}
}

void nb_cvfa_get_adj_graph(const nb_mesh2D_t *intmsh,
			   const nb_graph_t *trg_x_vol,
			   nb_graph_t *graph)
{
	graph->N = trg_x_vol->N;
	graph->wi = NULL;
	graph->wij = NULL;
	adj_graph_allocate_adj(graph, trg_x_vol, intmsh);
	adj_graph_set_adj(graph, trg_x_vol, intmsh);
}

static void adj_graph_allocate_adj(nb_graph_t *graph,
				   const nb_graph_t *trg_x_vol,
				   const nb_mesh2D_t *intmsh)
{
	uint32_t memsize_N_adj = graph->N * sizeof(*(graph->N_adj));
	uint32_t N_adj = adj_graph_get_N_adj(trg_x_vol, intmsh);
	uint32_t memsize_adj = graph->N * sizeof(*(graph->adj)) +
		N_adj * sizeof(**(graph->adj));
	char *memblock = nb_allocate_mem(memsize_N_adj + memsize_adj);
	graph->N_adj = (void*) memblock;
	graph->adj = (void*) (memblock + memsize_N_adj);	
}

static uint32_t adj_graph_get_N_adj(const nb_graph_t *trg_x_vol,
				    const nb_mesh2D_t *intmsh)
{
	nb_container_type cnt_type = NB_SORTED;
	uint32_t bank_size = nb_membank_get_memsize();
	uint32_t memsize = bank_size + nb_container_get_memsize(cnt_type);
	char *memblock = nb_soft_allocate_mem(memsize);
	nb_membank_t *membank = (void*) memblock;
	nb_container_t *list = (void*) (memblock + bank_size);

	nb_membank_init(membank, sizeof(uint32_t));

	nb_container_init(list, cnt_type);
	nb_container_set_comparer(list, compare_ids);

	uint32_t N = 0;
	for (uint32_t i = 0; i < trg_x_vol->N; i++) {
		adj_graph_get_list_x_vol(trg_x_vol, intmsh, i, membank, list);

		N += nb_container_get_length(list);

		while (nb_container_is_not_empty(list)) {
			uint32_t *id = nb_container_delete_first(list);
			nb_membank_free_mem(membank, id);
		}
	}
	
	nb_container_finish(list);
	nb_membank_finish(membank);
	nb_soft_free_mem(memsize, memblock);

	return N;
}

static void adj_graph_get_list_x_vol(const nb_graph_t *trg_x_vol,
				     const nb_mesh2D_t *intmsh,
				     uint32_t vol_id, nb_membank_t *bank,
				     nb_container_t *list)
{
	uint32_t N_vol_adj = trg_x_vol->N_adj[vol_id];
	for (uint32_t i = 0; i < N_vol_adj; i++) {
		uint32_t id = trg_x_vol->adj[vol_id][i];
		uint16_t N_adj = nb_mesh2D_elem_get_N_adj(intmsh, id);
		for (uint16_t j = 0; j < N_adj; j++) {
			uint32_t nid = nb_mesh2D_elem_get_adj(intmsh, id, j);
			if (vol_id != nid) {
				if (NULL == nb_container_exist(list, &nid)) {
					uint32_t *aux =
						nb_membank_allocate_mem(bank);
					*aux = nid;
					nb_container_insert(list, aux);
				}
			}
		}
	}
}

static void adj_graph_set_adj(nb_graph_t *graph,
			      const nb_graph_t *trg_x_vol,
			      const nb_mesh2D_t *intmsh)
{
	nb_container_type cnt_type = NB_SORTED;
	uint32_t bank_size = nb_membank_get_memsize();
	uint32_t memsize = bank_size + nb_container_get_memsize(cnt_type);
	char *memblock = nb_soft_allocate_mem(memsize);
	nb_membank_t *membank = (void*) memblock;
	nb_container_t *list = (void*) (memblock + bank_size);

	nb_membank_init(membank, sizeof(uint32_t));

	nb_container_init(list, cnt_type);
	nb_container_set_comparer(list, compare_ids);


	uint32_t mem_used = graph->N * sizeof(*(graph->N_adj)) +
		graph->N * sizeof(*(graph->adj));
	char *block = ((char*)graph->N_adj) + mem_used;

	for (uint32_t i = 0; i < trg_x_vol->N; i++) {
		adj_graph_get_list_x_vol(trg_x_vol, intmsh, i,
					  membank, list);

		uint32_t N_adj = nb_container_get_length(list);
		graph->N_adj[i] = N_adj;

		graph->adj[i] = (void*) block;
		block += N_adj * sizeof(**(graph->adj));

		uint16_t j = 0;
		while (nb_container_is_not_empty(list)) {
			uint32_t *id = nb_container_delete_first(list);

			graph->adj[i][j] = *id;
			j += 1;

			nb_membank_free_mem(membank, id);
		}
	}
	
	nb_container_finish(list);
	nb_membank_finish(membank);
	nb_soft_free_mem(memsize, memblock);
}

void nb_cvfa_load_trg_points(const nb_mesh2D_t *intmsh,
			     uint32_t trg_id, double t1[2],
			     double t2[2], double t3[2])
{
	uint32_t id1 = nb_mesh2D_elem_get_adj(intmsh, trg_id, 0);
	uint32_t id2 = nb_mesh2D_elem_get_adj(intmsh, trg_id, 1);
	uint32_t id3 = nb_mesh2D_elem_get_adj(intmsh, trg_id, 2);

	t1[0] = nb_mesh2D_node_get_x(intmsh, id1);
	t1[1] = nb_mesh2D_node_get_y(intmsh, id1);

	t2[0] = nb_mesh2D_node_get_x(intmsh, id2);
	t2[1] = nb_mesh2D_node_get_y(intmsh, id2);

	t3[0] = nb_mesh2D_node_get_x(intmsh, id3);
	t3[1] = nb_mesh2D_node_get_y(intmsh, id3);
}

void nb_cvfa_init_global_matrix(nb_sparse_t **K, const nb_graph_t *trg_x_vol,
				const nb_mesh2D_t *intmsh, int dof)
{
	uint32_t memsize = nb_graph_get_memsize();
	nb_graph_t *graph = nb_soft_allocate_mem(memsize);

	nb_graph_init(graph);
	nb_cvfa_get_adj_graph(intmsh, trg_x_vol, graph);

	/* Force symmetry for LU decomposition */
	nb_graph_force_symmetry(graph);

	*K = nb_sparse_create(graph, NULL, dof);

	nb_graph_finish(graph);
	nb_soft_free_mem(memsize, graph);
}

void nb_cvfa_load_faces(const nb_mesh2D_t *mesh,
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

void nb_cvfa_finish_faces(uint32_t N_faces, face_t **faces)
{
	for (uint32_t i = 0; i < N_faces; i++)
		nb_free_mem(faces[i]->subfaces);
}

bool nb_cvfa_face_is_internal(const face_t *face, const nb_mesh2D_t *mesh)
{
	uint32_t N_elems = nb_mesh2D_get_N_elems(mesh);
	return face->elems[1] < N_elems;
}

bool nb_cvfa_subface_in_simplex(const subface_t *subface)
{
	return (subface->N_int > 0);
}
