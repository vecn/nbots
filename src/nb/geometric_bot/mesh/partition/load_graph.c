#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

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
static uint16_t get_N_elems_around_2n(const nb_partition_t *part,
				      uint32_t elem_id, uint16_t face_id,
				      bool *is_interior_node);
static uint32_t get_right_elem(const nb_partition_t *part, uint32_t elem_id,
			       uint32_t ngb_id);
static void elem_by_nod_graph_count_adj(nb_graph_t *graph,
					const nb_partition_t *part);
static void elem_by_nod_graph_set_adj(nb_graph_t *graph,
				      const nb_partition_t *part);
static void set_elem_by_nod_adj(nb_graph_t *graph,
				const nb_partition_t *part,
				uint32_t elem_id);
static uint32_t set_elem_by_nod_rounding_adj(nb_graph_t *graph,
					     const nb_partition_t *part,
					     uint32_t elem_id,
					     uint16_t face_id,
					     uint16_t cnt);

void nb_partition_load_graph(const nb_partition_t *part,
			     vcn_graph_t *graph,
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
	char *memblock = malloc(memsize_N_adj + memsize_adj);
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
	char *memblock = malloc(memsize_N_adj + memsize_adj);
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
	char *memblock = malloc(memsize_N_adj + memsize_adj);
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
	char *memblock = malloc(memsize_N_adj + memsize_adj);
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
		if (nb_partition_elem_has_ngb(part, elem_id, j)) {
			bool is_interior_node;
			uint16_t N = get_N_elems_around_2n(part, elem_id, j,
							   &is_interior_node);
			if (is_interior_node)
				N_elem_adj += N - 2;
			else
				N_elem_adj += 1;
		}
	}
	return N_elem_adj;
}

static uint16_t get_N_elems_around_2n(const nb_partition_t *part,
				      uint32_t elem_id, uint16_t face_id,
				      bool *is_interior_node)
{
	uint32_t N_elems = nb_partition_get_N_elems(part);
	uint32_t N_adj = nb_partition_elem_get_N_adj(part, elem_id);

	uint16_t N = 0;
	uint32_t ngb = nb_partition_elem_get_ngb(part, elem_id,
						 (face_id + 1) % N_adj);
	if (ngb >= N_elems) {
		*is_interior_node = false;
		goto EXIT;
	}

	N = 1;
	uint32_t prev_id = ngb;
	uint32_t id = elem_id;
	while (ngb != id) {
		uint32_t aux = prev_id;
		prev_id = id;
		id = get_right_elem(part, id, aux);
		if (id >= N_elems) {
			*is_interior_node = false;
			goto EXIT;
		}
		N += 1;
	}
	*is_interior_node = true;
EXIT:
	return N;
}

static uint32_t get_right_elem(const nb_partition_t *part, uint32_t elem_id,
			       uint32_t ngb_id)
{
	uint32_t N_elems = nb_partition_get_N_elems(part);
	uint16_t N_adj = nb_partition_elem_get_N_adj(part, elem_id);
	uint16_t face_id = N_adj;
	for (uint16_t i = 0; i < N_adj; i++) {
		uint32_t ingb = nb_partition_elem_get_ngb(part, elem_id, i);
		if (ingb == ngb_id) {
			face_id = i;
			break;
		}
	}
	uint32_t right_elem;
	if (face_id < N_adj) {
		if (0 == face_id)
			face_id = N_adj - 1;
		else
			face_id -= 1;
		right_elem = nb_partition_elem_get_ngb(part, elem_id,
						       face_id);
	} else {
		right_elem = N_elems;
	}
	return right_elem;	
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
		if (nb_partition_elem_has_ngb(part, elem_id, j))
			cnt = set_elem_by_nod_rounding_adj(graph, part,
							   elem_id,
							   j, cnt);
	}
}

static uint32_t set_elem_by_nod_rounding_adj(nb_graph_t *graph,
					     const nb_partition_t *part,
					     uint32_t elem_id,
					     uint16_t face_id,
					     uint16_t cnt)
{
	bool is_interior_node;
	get_N_elems_around_2n(part, elem_id, face_id,
			      &is_interior_node);

	if (is_interior_node) {
		uint32_t N_adj = nb_partition_elem_get_N_adj(part, elem_id);
		uint32_t ngb = nb_partition_elem_get_ngb(part, elem_id,
							 (face_id+1) % N_adj);
		uint32_t prev_id = ngb;
		uint32_t id = elem_id;
		while (ngb != id) {
			uint32_t aux = prev_id;
			prev_id = id;
			id = get_right_elem(part, id, aux);
			if (ngb != id) {
				graph->adj[elem_id][cnt] = id;
				cnt += 1;				
			}
		}
	} else {
		uint32_t ngb = nb_partition_elem_get_ngb(part, elem_id,
							 face_id);
		graph->adj[elem_id][cnt] = ngb;
		cnt += 1;
	}
	return cnt;
}
