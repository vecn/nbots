#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "nb/memory_bot.h"
#include "nb/graph_bot.h"
#include "nb/geometric_bot/mesh/partition.h"
#include "nb/geometric_bot/mesh/partition/info.h"

static void load_graph_nodes_by_edges(const void *part, nb_graph_t *graph);
static void nodal_graph_allocate_adj(nb_graph_t *graph,
				     const nb_partition_t *part);
static void nodal_graph_count_adj(nb_graph_t *graph,
				  const nb_partition_t *part);
static void graph_assign_mem_adj(nb_graph_t *graph);
static void nodal_graph_set_adj(nb_graph_t *graph,
				const nb_partition_t *part);

static void load_graph_elems_by_edges(const void *part, nb_graph_t *graph);
static void elemental_graph_allocate_adj(nb_graph_t *graph,
					 const nb_partition_t *part);
static uint32_t get_N_elemental_adj(const nb_partition_t *part);
static uint16_t get_N_real_ngb(const nb_partition_t *part,
			       uint32_t elem_id);
static void elemental_graph_count_adj(nb_graph_t *graph,
				      const nb_partition_t *part);
static void elemental_graph_set_adj(nb_graph_t *graph,
				    const nb_partition_t *part);

static void load_graph_nodes_by_elems(const void *part, nb_graph_t *graph);
static void nod_by_elem_graph_allocate_adj(nb_graph_t *graph,
					   const nb_partition_t *part);
static uint32_t nod_by_elem_get_N_adj(const nb_partition_t *part);
static void nod_by_elem_graph_count_adj(nb_graph_t *graph,
					const nb_partition_t *part);
static void nod_by_elem_graph_set_adj(nb_graph_t *graph,
				      const nb_partition_t *part);
static void set_nod_by_elem_adj(nb_graph_t *graph, const nb_partition_t *part,
				uint32_t elem_id);

static void load_graph_elems_by_nodes(const void *part, nb_graph_t *graph);
static void elem_by_nod_graph_allocate_adj(nb_graph_t *graph,
					   const nb_partition_t *part);
static uint32_t elem_by_nod_get_N_adj(const nb_partition_t *part);
static uint32_t count_elem_by_nod_adj(const nb_partition_t *part,
				      uint32_t elem_id);
static uint16_t get_N_ngb_around_2n(const nb_partition_t *part,
				    uint32_t elem_id, uint16_t face_id);
static void elem_by_nod_graph_count_adj(nb_graph_t *graph,
					const nb_partition_t *part);
static void elem_by_nod_graph_set_adj(nb_graph_t *graph,
				      const nb_partition_t *part);
static void set_elem_by_nod_adj(nb_graph_t *graph,
				const nb_partition_t *part,
				uint32_t elem_id);
static uint16_t set_elem_by_nod_rounding_adj(nb_graph_t *graph,
					     const nb_partition_t *part,
					     uint32_t elem_id,
					     uint16_t face_id,
					     uint16_t cnt);
static void load_graph_elems_connected_to_nodes(const nb_partition_t *part,
						nb_graph_t *graph);
static void elem_to_nod_graph_allocate_adj(nb_graph_t *graph,
					   const nb_partition_t *part);
static uint32_t elem_to_nod_get_N_adj(const nb_partition_t *part);
static void elem_to_nod_graph_count_adj(nb_graph_t *graph,
					const nb_partition_t *part);
static void elem_to_nod_graph_set_adj(nb_graph_t *graph,
				      const nb_partition_t *part);

void nb_partition_load_graph(const nb_partition_t *part,
			     nb_graph_t *graph,
			     nb_partition_graph_type type)
{
	switch (type) {
	case NB_NODES_LINKED_BY_EDGES:
		load_graph_nodes_by_edges(part, graph);
		break;
	case NB_ELEMS_LINKED_BY_EDGES:
		load_graph_elems_by_edges(part, graph);
		break;
	case NB_NODES_LINKED_BY_ELEMS:
		load_graph_nodes_by_elems(part, graph);
		break;
	case NB_ELEMS_LINKED_BY_NODES:
		load_graph_elems_by_nodes(part, graph);
		break;
	case NB_ELEMS_CONNECTED_TO_NODES:
		load_graph_elems_connected_to_nodes(part, graph);
	default:
		graph->N = 0;
		break;
	}
}

static void load_graph_nodes_by_edges(const void *part, nb_graph_t *graph)
{
	graph->N = nb_partition_get_N_nodes(part);
	graph->wi = NULL;
	graph->wij = NULL;
	nodal_graph_allocate_adj(graph, part);
	nodal_graph_count_adj(graph, part);
	graph_assign_mem_adj(graph);
	nodal_graph_set_adj(graph, part);
}

static void nodal_graph_allocate_adj(nb_graph_t *graph,
				     const nb_partition_t *part)
{
	uint32_t N_edges = nb_partition_get_N_edges(part);

	uint32_t memsize_N_adj = graph->N * sizeof(*(graph->N_adj));
	uint32_t memsize_adj = 2 * N_edges * sizeof(**(graph->adj)) +
		graph->N * sizeof(*(graph->adj));
	char *memblock = nb_allocate_mem(memsize_N_adj + memsize_adj);
	graph->N_adj = (void*) memblock;
	graph->adj = (void*) (memblock + memsize_N_adj);
}

static void nodal_graph_count_adj(nb_graph_t *graph,
				  const nb_partition_t *part)
{
	memset(graph->N_adj, 0, graph->N * sizeof(*(graph->N_adj)));
	uint32_t N_edges = nb_partition_get_N_edges(part);
	for (uint32_t i = 0; i < N_edges; i++) {
		uint32_t id1 = nb_partition_edge_get_1n(part, i);
		uint32_t id2 = nb_partition_edge_get_2n(part, i);
		graph->N_adj[id1] += 1;
		graph->N_adj[id2] += 1;
	}
}

static void graph_assign_mem_adj(nb_graph_t *graph)
{
	uint32_t mem_used = graph->N * sizeof(*(graph->N_adj)) +
		graph->N * sizeof(*(graph->adj));
	char *block = (char*) graph->N_adj + mem_used;
	for (uint32_t i = 0; i < graph->N; i++) {
		graph->adj[i] = (void*) block;
		block += graph->N_adj[i] * sizeof(**(graph->adj));
	}
}

static void nodal_graph_set_adj(nb_graph_t *graph,
				const nb_partition_t *part)
{
	memset(graph->N_adj, 0, graph->N * sizeof(*(graph->N_adj)));
	uint32_t N_edges = nb_partition_get_N_edges(part);
	for (uint32_t i = 0; i < N_edges; i++) {
		uint32_t id1 = nb_partition_edge_get_1n(part, i);
		uint32_t id2 = nb_partition_edge_get_2n(part, i);
		graph->adj[id1][graph->N_adj[id1]] = id2;
		graph->adj[id2][graph->N_adj[id2]] = id1;
		graph->N_adj[id1] += 1;
		graph->N_adj[id2] += 1;
	}
}

static void load_graph_elems_by_edges(const void *part, nb_graph_t *graph)
{
	graph->N = nb_partition_get_N_elems(part);
	graph->wi = NULL;
	graph->wij = NULL;
	elemental_graph_allocate_adj(graph, part);
	elemental_graph_count_adj(graph, part);
	graph_assign_mem_adj(graph);
	elemental_graph_set_adj(graph, part);
}

static void elemental_graph_allocate_adj(nb_graph_t *graph,
					 const nb_partition_t *part)
{
	uint32_t memsize_N_adj = graph->N * sizeof(*(graph->N_adj));
	uint32_t N_adj = get_N_elemental_adj(part);
	uint32_t memsize_adj = graph->N * sizeof(*(graph->adj)) +
		N_adj * sizeof(**(graph->adj));
	char *memblock = nb_allocate_mem(memsize_N_adj + memsize_adj);
	graph->N_adj = (void*) memblock;
	graph->adj = (void*) (memblock + memsize_N_adj);	
}

static uint32_t get_N_elemental_adj(const nb_partition_t *part)
{
	uint32_t N_elems = nb_partition_get_N_elems(part);
	uint32_t N = 0;
	for (uint32_t i = 0; i < N_elems; i++)
		N += get_N_real_ngb(part, i);
	return N;
}

static uint16_t get_N_real_ngb(const nb_partition_t *part,
			       uint32_t elem_id)
{
	uint32_t N_adj = nb_partition_elem_get_N_adj(part, elem_id);
	uint16_t N_ngb = 0;
	for (uint32_t i = 0; i < N_adj; i++) {
		if (nb_partition_elem_has_ngb(part, elem_id, i))
			N_ngb += 1;
	}
	return N_ngb;
}

static void elemental_graph_count_adj(nb_graph_t *graph,
				      const nb_partition_t *part)
{
	uint32_t N_elems = nb_partition_get_N_elems(part);
	for (uint32_t i = 0; i < N_elems; i++)
		graph->N_adj[i] = get_N_real_ngb(part, i);
}

static void elemental_graph_set_adj(nb_graph_t *graph,
				    const nb_partition_t *part)
{
	uint32_t N_elems = nb_partition_get_N_elems(part);
	for (uint32_t i = 0; i < N_elems; i++) {
		uint16_t N_adj = nb_partition_elem_get_N_adj(part, i);
		uint16_t cnt = 0;
		for (uint32_t j = 0; j < N_adj; j++) {
			uint32_t ngb = nb_partition_elem_get_ngb(part, i, j);
			if (N_elems > ngb) {
				graph->adj[i][cnt] = ngb;
				cnt += 1;
			}
		}
	}
}

static void load_graph_nodes_by_elems(const void *part, nb_graph_t *graph)
{
	graph->N = nb_partition_get_N_nodes(part);
	graph->wi = NULL;
	graph->wij = NULL;
	nod_by_elem_graph_allocate_adj(graph, part);
	nod_by_elem_graph_count_adj(graph, part);
	graph_assign_mem_adj(graph);
	nod_by_elem_graph_set_adj(graph, part);
}

static void nod_by_elem_graph_allocate_adj(nb_graph_t *graph,
					   const nb_partition_t *part)
{
	uint32_t memsize_N_adj = graph->N * sizeof(*(graph->N_adj));
	uint32_t N_adj = nod_by_elem_get_N_adj(part);
	uint32_t memsize_adj = N_adj * sizeof(**(graph->adj)) +
		graph->N * sizeof(*(graph->adj));
	char *memblock = nb_allocate_mem(memsize_N_adj + memsize_adj);
	graph->N_adj = (void*) memblock;
	graph->adj = (void*) (memblock + memsize_N_adj);
}

static uint32_t nod_by_elem_get_N_adj(const nb_partition_t *part)
{
	uint32_t N_edges = nb_partition_get_N_edges(part);
	uint32_t N = 2 * N_edges;

	uint32_t N_elems = nb_partition_get_N_elems(part);
	for (uint32_t i = 0; i < N_elems; i++) {
		uint16_t N_adj = nb_partition_elem_get_N_adj(part, i);
		N += N_adj * (N_adj - 3);
	}
	return N;
}

static void nod_by_elem_graph_count_adj(nb_graph_t *graph,
					const nb_partition_t *part)
{
	uint32_t N_edges = nb_partition_get_N_edges(part);
	memset(graph->N_adj, 0, graph->N * sizeof(*(graph->N_adj)));
	for (uint32_t i = 0; i < N_edges; i++) {
		uint32_t id1 = nb_partition_edge_get_1n(part, i);
		uint32_t id2 = nb_partition_edge_get_2n(part, i);
		graph->N_adj[id1] += 1;
		graph->N_adj[id2] += 1;
	}

	uint32_t N_elems = nb_partition_get_N_elems(part);
	for (uint32_t i = 0; i < N_elems; i++) {
		uint16_t N_adj = nb_partition_elem_get_N_adj(part, i);
		if (3 < N_adj) {
			for (uint32_t j = 0; j < N_adj; j++) {
				uint32_t id = 
					nb_partition_elem_get_adj(part,	i, j);
				graph->N_adj[id] += N_adj - 3;
			}
		}
	}
}

static void nod_by_elem_graph_set_adj(nb_graph_t *graph,
				      const nb_partition_t *part)
{
	uint32_t N_edges = nb_partition_get_N_edges(part);
	memset(graph->N_adj, 0, graph->N * sizeof(*(graph->N_adj)));
	for (uint32_t i = 0; i < N_edges; i++) {
		uint32_t id1 = nb_partition_edge_get_1n(part, i);
		uint32_t id2 = nb_partition_edge_get_2n(part, i);
		graph->adj[id1][graph->N_adj[id1]] = id2;
		graph->adj[id2][graph->N_adj[id2]] = id1;
		graph->N_adj[id1] += 1;
		graph->N_adj[id2] += 1;
	}
	uint32_t N_elems = nb_partition_get_N_elems(part);
	for (uint32_t i = 0; i < N_elems; i++)
		set_nod_by_elem_adj(graph, part, i);
}

static void set_nod_by_elem_adj(nb_graph_t *graph, const nb_partition_t *part,
				uint32_t elem_id)
{
	uint16_t N_adj = nb_partition_elem_get_N_adj(part, elem_id);
	if (3 >= N_adj)
		goto EXIT;

	for (uint16_t i = 0; i < N_adj; i++) {
		uint32_t id1 = nb_partition_elem_get_adj(part, elem_id, i);
		for (uint16_t j = 0; j < N_adj - 3; j++) {
			uint32_t k = (i + j + 2) % N_adj;
			uint32_t id2 = 
				nb_partition_elem_get_adj(part, elem_id, k);
			graph->adj[id1][graph->N_adj[id1]] = id2;
			graph->N_adj[id1] += 1;
		}
	}
EXIT:
	return;
}

static void load_graph_elems_by_nodes(const void *part, nb_graph_t *graph)
{
	graph->N = nb_partition_get_N_elems(part);
	graph->wi = NULL;
	graph->wij = NULL;
	elem_by_nod_graph_allocate_adj(graph, part);
	elem_by_nod_graph_count_adj(graph, part);
	graph_assign_mem_adj(graph);
	elem_by_nod_graph_set_adj(graph, part);
}

static void elem_by_nod_graph_allocate_adj(nb_graph_t *graph,
					   const nb_partition_t *part)
{
	uint32_t memsize_N_adj = graph->N * sizeof(*(graph->N_adj));
	uint32_t N_adj = elem_by_nod_get_N_adj(part);
	uint32_t memsize_adj = graph->N * sizeof(*(graph->adj)) +
		N_adj * sizeof(**(graph->adj));
	char *memblock = nb_allocate_mem(memsize_N_adj + memsize_adj);
	graph->N_adj = (void*) memblock;
	graph->adj = (void*) (memblock + memsize_N_adj);	
}

static uint32_t elem_by_nod_get_N_adj(const nb_partition_t *part)
{
	uint32_t N = 0;
	uint32_t N_elems = nb_partition_get_N_elems(part);
	for (uint32_t i = 0; i < N_elems; i++)
		N += count_elem_by_nod_adj(part, i);
	return N;
}

static uint32_t count_elem_by_nod_adj(const nb_partition_t *part,
				      uint32_t elem_id)
{
	uint32_t N_elem_adj = 0;
	uint16_t N_adj = nb_partition_elem_get_N_adj(part, elem_id);
	for (uint16_t j = 0; j < N_adj; j++) {
		uint16_t N = get_N_ngb_around_2n(part, elem_id, j);
		N_elem_adj += N;
	}
	return N_elem_adj;
}

static uint16_t get_N_ngb_around_2n(const nb_partition_t *part,
				    uint32_t elem_id, uint16_t face_id)
{
	uint32_t N_elems = nb_partition_get_N_elems(part);

	uint16_t N = 0;

	uint32_t front_ngb_id = nb_partition_elem_get_ngb(part, elem_id, face_id);
	uint32_t nid_prev = elem_id;
	uint32_t nid = nb_partition_elem_face_get_right_ngb(part, elem_id,
							    face_id);
	while(nid != front_ngb_id && nid < N_elems) {
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
		nid = nb_partition_elem_face_get_left_ngb(part, front_ngb_id, aux);
		while(nid < N_elems) {
			N += 1;
			uint16_t aux = nb_partition_elem_ngb_get_face(part, nid,
								      nid_prev);
			nid_prev = nid;
			nid = nb_partition_elem_face_get_left_ngb(part, nid,
								  aux);
		}

	}
	return N;
}

static void elem_by_nod_graph_count_adj(nb_graph_t *graph,
					const nb_partition_t *part)
{
	uint32_t N_elems = nb_partition_get_N_elems(part);
	for (uint32_t i = 0; i < N_elems; i++)
		graph->N_adj[i] = count_elem_by_nod_adj(part, i);
}

static void elem_by_nod_graph_set_adj(nb_graph_t *graph,
				      const nb_partition_t *part)
{
	uint32_t N_elems = nb_partition_get_N_elems(part);
	for (uint32_t i = 0; i < N_elems; i++)
		set_elem_by_nod_adj(graph, part, i);
}

static void set_elem_by_nod_adj(nb_graph_t *graph,
				const nb_partition_t *part,
				uint32_t elem_id)
{
	uint16_t cnt = 0;
	uint16_t N_adj = nb_partition_elem_get_N_adj(part, elem_id);
	for (uint16_t j = 0; j < N_adj; j++) {
		cnt = set_elem_by_nod_rounding_adj(graph, part,
						   elem_id,
						   j, cnt);
	}
}

static uint16_t set_elem_by_nod_rounding_adj(nb_graph_t *graph,
					     const nb_partition_t *part,
					     uint32_t elem_id,
					     uint16_t face_id,
					     uint16_t cnt)
{
	uint32_t N_elems = nb_partition_get_N_elems(part);
	uint16_t N_adj = nb_partition_elem_get_N_adj(part, elem_id);

	uint32_t nid_prev = elem_id;
	uint32_t nid = nb_partition_elem_get_ngb(part, elem_id, face_id);
	uint32_t front_ngb = nb_partition_elem_get_ngb(part, elem_id,
						       (face_id+1)%N_adj);
	while(nid != front_ngb && nid < N_elems) {
		graph->adj[elem_id][cnt] = nid;
		cnt += 1;
		uint16_t aux = nb_partition_elem_ngb_get_face(part, nid,
							      nid_prev);
		nid_prev = nid;
		nid = nb_partition_elem_face_get_right_ngb(part, nid, aux);
	}
	if (nid >= N_elems && front_ngb < N_elems) {
		uint16_t aux = nb_partition_elem_ngb_get_face(part, front_ngb,
							      elem_id);
		uint32_t nid_prev = front_ngb;
		uint32_t nid = nb_partition_elem_face_get_left_ngb(part,
								   front_ngb,
								   aux);
		while(nid < N_elems) {
			graph->adj[elem_id][cnt] = nid;
			cnt += 1;
			aux = nb_partition_elem_ngb_get_face(part, nid,
							     nid_prev);
			nid_prev = nid;
			nid = nb_partition_elem_face_get_left_ngb(part, nid,
								  aux);
		}

	}
	return cnt;
}

static void load_graph_elems_connected_to_nodes(const nb_partition_t *part,
						nb_graph_t *graph)
{
	graph->N = nb_partition_get_N_nodes(part);
	graph->wi = NULL;
	graph->wij = NULL;
	elem_to_nod_graph_allocate_adj(graph, part);
	elem_to_nod_graph_count_adj(graph, part);
	graph_assign_mem_adj(graph);
	elem_to_nod_graph_set_adj(graph, part);
}

static void elem_to_nod_graph_allocate_adj(nb_graph_t *graph,
					   const nb_partition_t *part)
{
	uint32_t memsize_N_adj = graph->N * sizeof(*(graph->N_adj));
	uint32_t N_adj = elem_to_nod_get_N_adj(part);
	uint32_t memsize_adj = graph->N * sizeof(*(graph->adj)) +
		N_adj * sizeof(**(graph->adj));
	char *memblock = nb_allocate_mem(memsize_N_adj + memsize_adj);
	graph->N_adj = (void*) memblock;
	graph->adj = (void*) (memblock + memsize_N_adj);	
}

static uint32_t elem_to_nod_get_N_adj(const nb_partition_t *part)
{
	uint32_t N_elems = nb_partition_get_N_elems(part);
	uint32_t N;
	if (NB_TRIAN == nb_partition_get_type(part)) {
		N = 3 * N_elems;
	} else {
		N = 0;
		for (uint32_t i = 0; i < N_elems; i++)
			N += nb_partition_elem_get_N_adj(part, i);
	}
	return N;
}

static void elem_to_nod_graph_count_adj(nb_graph_t *graph,
					const nb_partition_t *part)
{
	uint32_t N_elems = nb_partition_get_N_elems(part);

	memset(graph->N_adj, 0, graph->N * sizeof(*(graph->N_adj)));
	for (uint32_t i = 0; i < N_elems; i++) {
		uint16_t N_adj = nb_partition_elem_get_N_adj(part, i);
		for (uint16_t j = 0; j < N_adj; j++) {
			uint32_t id = nb_partition_elem_get_adj(part, i, j);
			graph->N_adj[id] += 1;
		}
	}
}

static void elem_to_nod_graph_set_adj(nb_graph_t *graph,
				      const nb_partition_t *part)
{
	uint32_t N_elems = nb_partition_get_N_elems(part);

	memset(graph->N_adj, 0, graph->N * sizeof(*(graph->N_adj)));
	for (uint32_t i = 0; i < N_elems; i++) {
		uint16_t N_adj = nb_partition_elem_get_N_adj(part, i);
		for (uint16_t j = 0; j < N_adj; j++) {
			uint32_t id = nb_partition_elem_get_adj(part, i, j);
			uint32_t cnt = graph->N_adj[id];
			graph->adj[id][cnt] = i;
			graph->N_adj[id] = cnt + 1;
		}
	}
}
