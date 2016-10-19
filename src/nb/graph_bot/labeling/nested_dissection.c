#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>

#include "nb/memory_bot.h"
#include "nb/graph_bot/graph.h"
#include "nb/graph_bot/labeling/nested_dissection.h"

static uint32_t get_N_total_adj(const nb_graph_t *graph);

void nb_graph_labeling_nd(const nb_graph_t *graph,
			  uint32_t *perm, uint32_t *iperm)
{
    ;
}

static uint32_t get_N_total_adj(const nb_graph_t *graph)
{
	uint32_t N = 0;
	for (uint32_t i = 0; i < graph->N; i++)
		N += graph->N_adj[i];
	return N;
}

