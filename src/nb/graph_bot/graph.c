#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>

#include "nb/memory_bot.h"
#include "nb/graph_bot/graph.h"
#include "nb/graph_bot/labeling/amd.h"


static void allocate_and_copy_graph(nb_graph_t *graph,
				    const nb_graph_t *graph_src);
static uint32_t get_N_total_adj(const nb_graph_t *graph);
static void copy_edge_weights(nb_graph_t *graph, const nb_graph_t *graph_src);
static uint32_t get_wij_memsize(nb_graph_t *graph);
static void graph_clear(nb_graph_t *graph);

uint32_t nb_graph_get_memsize(void)
{
	return sizeof(nb_graph_t);
}

void nb_graph_init(nb_graph_t *graph)
{
	memset(graph, 0, nb_graph_get_memsize());
}

void nb_graph_copy(void *graph_ptr, const void *graph_src_ptr)
{
	nb_graph_t *graph = graph_ptr;
	const nb_graph_t *graph_src = graph_src_ptr;

	memset(graph, 0, nb_graph_get_memsize());
	graph->N = graph_src->N;

	if (NULL != graph_src->N_adj)
		allocate_and_copy_graph(graph, graph_src);

	if (NULL != graph_src->wi) {
		nb_graph_init_vtx_weights(graph);
		memcpy(graph->wi, graph_src->wi,
		       graph->N * sizeof(*(graph->wi)));
	}

	if (NULL != graph_src->wij) {
		nb_graph_init_edge_weights(graph);
		copy_edge_weights(graph, graph_src);
	}
}

static void allocate_and_copy_graph(nb_graph_t *graph,
				    const nb_graph_t *graph_src)
{
	uint32_t N = graph_src->N;
	uint32_t N_total_adj = get_N_total_adj(graph_src);
	uint32_t memsize = N * (sizeof(*(graph->N_adj)) +
				sizeof(*(graph->adj))) +
		N_total_adj * sizeof(**(graph->adj));
	char *memblock = nb_allocate_mem(memsize);
	memset(memblock, 0, memsize);

	graph->N_adj = (void*) memblock;
	graph->adj = (void*) (memblock + N * sizeof(*(graph->N_adj)));

	memblock += N * (sizeof(*(graph->N_adj)) +
			 sizeof(*(graph->adj)));

	for (uint32_t i = 0; i < graph_src->N; i++) {
		uint32_t N_adj = graph_src->N_adj[i];
		graph->N_adj[i] = N_adj;
		graph->adj[i] = (void*) memblock;
		uint32_t memstep = N_adj * sizeof(**(graph->adj));
		memblock += memstep;
		memcpy(graph->adj[i], graph_src->adj[i], memstep);
	}
}

static uint32_t get_N_total_adj(const nb_graph_t *graph)
{
	uint32_t N = 0;
	for (uint32_t i = 0; i < graph->N; i++)
		N += graph->N_adj[i];
	return N;
}

static void copy_edge_weights(nb_graph_t *graph, const nb_graph_t *graph_src)
{
	for (uint32_t i = 0; i < graph->N; i++)
		memcpy(graph->wij[i], graph_src->wij[i],
		       graph->N_adj[i] * sizeof(**(graph->wij)));
}

void nb_graph_finish(nb_graph_t *graph)
{
	graph_clear(graph);
}

static void graph_clear(nb_graph_t *graph)
{
	if (NULL != graph->N_adj)
		/* Includes graph->adj in the same memblock */
		nb_free_mem(graph->N_adj);

	if (NULL != graph->wi)
		nb_free_mem(graph->wi);

	if (NULL != graph->wij)
		/* Includes graph->wij[i] in the same memblock */
		nb_free_mem(graph->wij);
}

void nb_graph_clear(nb_graph_t *graph)
{
	graph_clear(graph);
	memset(graph, 0, nb_graph_get_memsize());	
}

void nb_graph_init_vtx_weights(nb_graph_t *graph)
{
	graph->wi = nb_allocate_zero_mem(graph->N * sizeof(*(graph->wi)));
}

void nb_graph_init_edge_weights(nb_graph_t *graph)
{
	uint32_t memsize = get_wij_memsize(graph);
	char *memblock = nb_allocate_zero_mem( memsize);
	graph->wij = (void*) memblock;
	char *block =  memblock + graph->N * sizeof(*(graph->wij));
	for (uint32_t i = 0; i < graph->N; i++) {
		graph->wij[i] = (void*) block;
		block += graph->N_adj[i] * sizeof(**(graph->wij));
	}
}

static uint32_t get_wij_memsize(nb_graph_t *graph)
{
	uint32_t size = graph->N * sizeof(*(graph->wij));
	for (uint32_t i = 0; i < graph->N; i++)
		size += graph->N_adj[i] * sizeof(**(graph->wij));
	return size;
}

void nb_graph_finish_vtx_weights(nb_graph_t *graph)
{
	nb_free_mem(graph->wi);
	graph->wi = NULL;
}

void nb_graph_finish_edge_weights(nb_graph_t *graph)
{
	nb_free_mem(graph->wij);
	graph->wij = NULL;
}

uint32_t nb_graph_find_N_intersected_adj(const nb_graph_t *graph,
					 uint32_t id1, uint32_t id2)
{
	uint32_t N = 0;
	for (uint32_t i = 0; i < graph->N_adj[id1]; i++) {
		for (uint32_t j = 0; j < graph->N_adj[id2]; j++) {
			if (graph->adj[id1][i] == graph->adj[id2][j]) {
				N += 1;
				break;
			}
		}	
	}
	return N;
}

void nb_graph_labeling(const nb_graph_t *const graph,
		       uint32_t *perm, uint32_t* iperm,
		       nb_labeling_algorithm alg)
{
	switch (alg) {
	case NB_LABELING_MMD:
		nb_graph_labeling_mmd(graph, perm, iperm);
		break;
	case NB_LABELING_AMD:
		nb_graph_labeling_amd(graph, perm, iperm);
		break;
	case NB_LABELING_ND:
		nb_graph_labeling_nd(graph, perm, iperm);
		break;
	default:
		nb_graph_labeling_amd(graph, perm, iperm);
	}
}
