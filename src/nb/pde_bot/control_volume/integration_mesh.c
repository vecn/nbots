#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>

#include "nb/memory_bot.h"
#include "nb/container_bot.h"
#include "nb/geometric_bot.h"

#include "integration_mesh.h"

#define INTEGRATOR_TYPE NB_TRIAN

static void init_containers_trg_x_vol(nb_container_t **all_trg_x_vol,
				      nb_container_type cnt_type,
				      nb_membank_t *membank,
				      const nb_partition_t *part,
				      const nb_partition_t *intmsh);
static void vol_get_adj(const nb_partition_t *part,
			const nb_partition_t *intmsh,
			const nb_graph_t *trg_x_vtx,
			nb_membank_t *membank,
			uint32_t vol_id,
			nb_container_t *trg_adj);
static int8_t compare_ids(const void *ptr1, const void *ptr2);
static bool vol_intersects_trg(const nb_partition_t *part,
			       const nb_partition_t *intmsh,
			       uint32_t vol_id, uint32_t trg_id);
static void mesh_load_sgm_from_adj(const nb_partition_t *part,
				   uint32_t elem_id, uint16_t adj_id,
				   double s1[2], double s2[2]);
static void put_neighbours_in_active(const nb_partition_t *intmsh,
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

uint32_t nb_cvfa_get_integration_mesh_memsize(void)
{
	return nb_partition_get_memsize(INTEGRATOR_TYPE);
}

void nb_cvfa_init_integration_mesh(nb_partition_t *intmsh)
{
	nb_partition_init(intmsh , INTEGRATOR_TYPE);
}

void nb_cvfa_load_integration_mesh(const nb_partition_t *part,
				   nb_partition_t *intmsh)
{
	uint32_t N_elems = nb_partition_get_N_elems(part);

	uint32_t mesh_size = nb_mesh_get_memsize();
	uint32_t vtx_size = 2 * N_elems * sizeof(double);
	uint32_t perm_size = N_elems * sizeof(uint32_t);
	uint32_t memsize = mesh_size + vtx_size + perm_size;
	char *memblock = NB_SOFT_MALLOC(memsize);
	nb_mesh_t *mesh = (void*) memblock;
	double *vtx = (void*) (memblock + mesh_size);
	uint32_t *perm = (void*) (memblock + mesh_size + vtx_size);
	
	for (uint32_t i = 0; i < N_elems; i++) {
		vtx[i * 2] = nb_partition_elem_get_x(part, i);
		vtx[i*2+1] = nb_partition_elem_get_y(part, i);
	}	

	nb_mesh_init(mesh);
	nb_mesh_get_smallest_ns_alpha_complex(mesh, N_elems, vtx, 0.666);
	nb_partition_load_from_mesh(intmsh, mesh);
	nb_mesh_finish(mesh);

	for (uint32_t i = 0; i < N_elems; i++) {
		uint32_t id = nb_partition_get_invtx(intmsh, i);
		perm[id] = i;
	}

	nb_partition_set_nodal_permutation(intmsh, perm);

	NB_SOFT_FREE(memsize, memblock);
}

void nb_cvfa_correlate_partition_and_integration_mesh
					(const nb_partition_t *part,
					 const nb_partition_t *intmsh,
					 nb_graph_t *trg_x_vol)
{
	nb_container_type cnt_type = NB_SORTED;
	uint32_t N_elems = nb_partition_get_N_elems(part);
	uint32_t cnt_size = nb_container_get_memsize(cnt_type);
	uint32_t bank_size = nb_membank_get_memsize();
	uint32_t memsize = N_elems * (cnt_size + sizeof(void*)) + bank_size;
	char *memblock = NB_SOFT_MALLOC(memsize);
	nb_container_t **all_trg_x_vol = (void*) memblock;
	nb_membank_t *membank = (void*) (memblock + N_elems * 
					 (cnt_size + sizeof(void*)));

	nb_membank_init(membank, sizeof(uint32_t));

	init_containers_trg_x_vol(all_trg_x_vol, cnt_type, membank,
				  part, intmsh);
	
	trg_x_vol->N = N_elems;
	trg_x_vol_allocate_adj(trg_x_vol, all_trg_x_vol);
	trg_x_vol_set_adj(trg_x_vol, all_trg_x_vol, membank);

	finish_containers_trg_x_vol(N_elems, all_trg_x_vol, membank);

	nb_membank_finish(membank);
	NB_SOFT_FREE(memsize, memblock);
}

static void init_containers_trg_x_vol(nb_container_t **all_trg_x_vol,
				      nb_container_type cnt_type,
				      nb_membank_t *membank,
				      const nb_partition_t *part,
				      const nb_partition_t *intmsh)
{
	uint32_t memsize = nb_graph_get_memsize();
	char *memblock = NB_SOFT_MALLOC(memsize);
	nb_graph_t *trg_x_vtx = (void*) memblock;

	nb_graph_init(trg_x_vtx);
	nb_partition_load_graph(part, trg_x_vtx, NB_ELEMS_CONNECTED_TO_NODES);

	uint32_t N = 0;
	uint32_t N_elems = nb_partition_get_N_elems(part);
	char *block = ((char*) all_trg_x_vol) + N_elems * sizeof(void*);
	uint32_t cnt_size = nb_container_get_memsize(cnt_type);
	for (uint32_t i = 0; i < N_elems; i++) {
		nb_container_t *trg_adj = (void*) block;

		all_trg_x_vol[i] = trg_adj;
		nb_container_init(trg_adj, cnt_type);

		printf(" -- FINDING BUG 1 %i\n", i);/* TEMPORAL AQUI VOY */
		vol_get_adj(part, intmsh, trg_x_vtx, membank, i, trg_adj);
		printf(" -- FINDING BUG 3\n");/* TEMPORAL */

		block += cnt_size;
	}

	nb_graph_finish(trg_x_vtx);

	NB_SOFT_FREE(memsize, memblock);
}

static void vol_get_adj(const nb_partition_t *part,
			const nb_partition_t *intmsh,
			const nb_graph_t *trg_x_vtx,
			nb_membank_t *membank,
			uint32_t vol_id,
			nb_container_t *trg_adj)
{
	uint32_t cnt_size = nb_container_get_memsize(NB_SORTED);
	uint32_t active_size = nb_container_get_memsize(NB_QUEUE);
	uint32_t memsize = cnt_size + active_size;
	char *memblock = NB_SOFT_MALLOC(memsize);
	nb_container_t *out = (void*) memblock;
	nb_container_t *active = (void*) (memblock + cnt_size);

	nb_container_init(out, NB_SORTED);
	nb_container_init(active, NB_QUEUE);

	nb_container_set_comparer(trg_adj, compare_ids);
	nb_container_set_comparer(out, compare_ids);

	uint32_t *id = nb_membank_allocate_mem(membank);
	*id = trg_x_vtx->adj[vol_id][0];
	nb_container_insert(active, id);
	
	while (nb_container_is_not_empty(active)) {
		uint32_t *id = nb_container_delete_first(active);
		bool intersection = vol_intersects_trg(part, intmsh,
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
		uint32_t *id = nb_container_delete_first(out);
		nb_membank_free_mem(membank, id);
	}

	nb_container_finish(out);
	nb_container_finish(active);

	NB_SOFT_FREE(memsize, memblock);
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

static bool vol_intersects_trg(const nb_partition_t *part,
			       const nb_partition_t *intmsh,
			       uint32_t vol_id, uint32_t trg_id)
{
	bool out = false;

	double a1[2], a2[2], b1[2], b2[2];
	
	uint16_t N_adj1 = nb_partition_elem_get_N_adj(part, vol_id);
	for (uint16_t i = 0; i < N_adj1; i++) {
		mesh_load_sgm_from_adj(part, vol_id, i, a1, a2);
		uint16_t N_adj2 = nb_partition_elem_get_N_adj(intmsh, trg_id);
		for (uint16_t j = 0; j < N_adj2; j++) {
			mesh_load_sgm_from_adj(intmsh, trg_id, j, b1, b2);
			out = vcn_utils2D_are_sgm_intersected(a1, a2, b1,
							      b2, NULL);
			if (out)
				goto EXIT;
		}
	}
EXIT:
	return out;
}

static void mesh_load_sgm_from_adj(const nb_partition_t *part,
				   uint32_t elem_id, uint16_t adj_id,
				   double s1[2], double s2[2])
{
	uint32_t N_adj = nb_partition_elem_get_N_adj(part, elem_id);
	uint32_t id1 = nb_partition_elem_get_adj(part, elem_id, adj_id);
	uint32_t id2 = nb_partition_elem_get_adj(part, elem_id,
						 (adj_id + 1) % N_adj);
	s1[0] = nb_partition_node_get_x(part, id1);
	s1[1] = nb_partition_node_get_y(part, id1);

	s2[0] = nb_partition_node_get_x(part, id2);
	s2[1] = nb_partition_node_get_y(part, id2);
}

static void put_neighbours_in_active(const nb_partition_t *intmsh,
				     nb_membank_t *membank,
				     const nb_container_t *trg_adj,
				     const nb_container_t *out,
				     nb_container_t *active, uint32_t id)
{
	uint16_t N_adj = nb_partition_elem_get_N_adj(intmsh, id);
	for (uint32_t i = 0; i < N_adj; i++) {
		if (nb_partition_elem_has_ngb(intmsh, id, i)) {
			uint32_t ngb = nb_partition_elem_get_ngb(intmsh,
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
			uint32_t *id = nb_membank_allocate_mem(membank);
			*id = ngb;
			nb_container_insert(active, id);
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
