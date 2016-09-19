#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>

#include "nb/graph_bot/graph.h"
#include "nb/graph_bot/labeling/amd.h"

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

void nb_graph_finish(nb_graph_t *graph)
{
	graph_clear(graph);
}

static void graph_clear(nb_graph_t *graph)
{
	if (NULL != graph->N_adj)
		/* Includes graph->adj in the same memblock */
		free(graph->N_adj);

	if (NULL != graph->wi)
		free(graph->wi);

	if (NULL != graph->wij)
		/* Includes graph->wij[i] in the same memblock */
		free(graph->wij);
}

void nb_graph_clear(nb_graph_t *graph)
{
	graph_clear(graph);
	memset(graph, 0, nb_graph_get_memsize());	
}

nb_graph_t* nb_graph_create(void)
{
	return calloc(1, sizeof(nb_graph_t));
}

void nb_graph_destroy(nb_graph_t *graph)
{
	nb_graph_finish(graph);
	free(graph);
}

void nb_graph_init_vtx_weights(nb_graph_t *graph)
{
	graph->wi = calloc(graph->N, sizeof(*(graph->wi)));
}

void nb_graph_init_edge_weights(nb_graph_t *graph)
{
	uint32_t memsize = get_wij_memsize(graph);
	char *memblock = calloc(1, memsize);
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
	free(graph->wi);
	graph->wi = NULL;
}

void nb_graph_finish_edge_weights(nb_graph_t *graph)
{
	free(graph->wij);
	graph->wij = NULL;
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
	default:
		nb_graph_labeling_amd(graph, perm, iperm);
	}
}
