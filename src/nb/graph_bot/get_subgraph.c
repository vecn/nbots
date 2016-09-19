#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>

#include "nb/graph_bot/graph.h"

static uint32_t exist(uint32_t N, uint32_t *array, uint32_t val);
static void graph_clear(nb_graph_t *graph);

static uint32_t exist(uint32_t N, uint32_t *array, uint32_t val)
{
	/* TEMPORAL Exahustive search */
	for (uint32_t i = 0; i < N; i++) {
		if (array[i] == val)
			return i;
	}
	return N;
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


nb_graph_t* nb_graph_get_subgraph(const nb_graph_t *const graph,
				    uint32_t N_nodes, uint32_t *nodes)
/* DEPRECATED: Check memory allocation */
{
	nb_graph_t *subgraph = calloc(1, sizeof(nb_graph_t));

	subgraph->N = N_nodes;
	uint32_t memsize = N_nodes * 
		(sizeof(*(subgraph->N_adj)) + sizeof(*(subgraph->adj)));
	char *memblock = malloc(memsize);
	subgraph->N_adj = (void*) memblock;
	subgraph->adj = (void*) (memblock + N_nodes *
				 sizeof(*(subgraph->N_adj)));
	memset(subgraph->adj, 0, N_nodes * sizeof(*(subgraph->adj)));

	if (NULL != graph->wi)
		subgraph->wi = malloc(N_nodes * sizeof(*(subgraph->wi)));

	if (NULL != graph->wij)
		subgraph->wij = calloc(N_nodes, sizeof(*(subgraph->wij)));

	for (uint32_t i = 0; i < N_nodes; i++) {
		subgraph->N_adj[i] = 0;
		for (uint32_t j = 0; j < graph->N_adj[nodes[i]]; j++) {
			uint32_t id = exist(N_nodes, nodes, graph->adj[nodes[i]][j]);
			if (id < N_nodes)
				subgraph->N_adj[i] += 1;
		}
		if (subgraph->N_adj[i] > 0) {
			subgraph->adj[i] =
				malloc(subgraph->N_adj[i] * sizeof(*(subgraph->adj[i])));
			if (NULL != graph->wij)
				subgraph->wij[i] = 
					malloc(subgraph->N_adj[i] * sizeof(*(subgraph->wij[i])));
			uint32_t cnt = 0;
			for (uint32_t j = 0; j < graph->N_adj[nodes[i]]; j++) {
				uint32_t id = exist(N_nodes, nodes, graph->adj[nodes[i]][j]);
				if (id < N_nodes) {
					subgraph->adj[i][cnt] = id;
					if (NULL != graph->wij)
						subgraph->wij[i][cnt] = graph->wij[nodes[i]][j];
					cnt += 1;
				}
			}
		}
		if (NULL != graph->wi) 
			subgraph->wi[i] = graph->wi[nodes[i]];
  }

  return subgraph;
}
