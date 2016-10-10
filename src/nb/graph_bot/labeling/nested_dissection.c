#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <alloca.h>

#include "../imported_libs/metis.h"

#include "nb/memory_bot.h"
#include "nb/graph_bot/graph.h"
#include "nb/graph_bot/labeling/nested_dissection.h"

static uint32_t get_N_total_adj(const nb_graph_t *graph);
static void get_crs(const nb_graph_t *graph, idx_t* xadj, idx_t *adj);

void nb_graph_labeling_nd(const nb_graph_t *graph,
			  uint32_t *perm, uint32_t *iperm)
{
	idx_t N = graph->N;

	idx_t op[METIS_NOPTIONS];
	METIS_SetDefaultOptions(op);
	op[METIS_OPTION_NUMBERING] = 0;

	uint32_t N_total_adj = get_N_total_adj(graph);

	uint32_t memsize = (1 + 3 * N + N_total_adj) * sizeof(idx_t);
	char *memblock = nb_soft_allocate_mem(memsize);
	
	idx_t *xadj = (void*) memblock;
	idx_t *adj = (void*) (memblock + (1 + N) * sizeof(idx_t));
	idx_t *metis_perm = (void*)
		(memblock + (1 + N + N_total_adj) * sizeof(idx_t));
	idx_t *metis_iperm = (void*)
		(memblock + (1 + 2 * N + N_total_adj) * sizeof(idx_t));

	get_crs(graph, xadj, adj);

	int status = METIS_NodeND(&N, xadj, adj, NULL, op,
				  metis_perm, metis_iperm);
	if (METIS_OK != status) {
		fprintf(stderr, "GRAPH_BOT: METIS_NodeND() fails\n");
		exit(1);
	}

	for (uint32_t i = 0; i < N; i++) {
		perm[i] = metis_perm[i];
		iperm[i] = metis_iperm[i];
	}

	nb_soft_free_mem(memsize, memblock);
}

static uint32_t get_N_total_adj(const nb_graph_t *graph)
{
	uint32_t N = 0;
	for (uint32_t i = 0; i < graph->N; i++)
		N += graph->N_adj[i];
	return N;
}

static void get_crs(const nb_graph_t *graph, idx_t* xadj, idx_t *adj)
{
	uint32_t id = 0;
	for (uint32_t i = 0; i < graph->N; i++) {
		xadj[i] = id;
		uint32_t N_adj = graph->N_adj[i];
		for (uint16_t j = 0; j < N_adj; j++) {
			adj[id] = graph->adj[i][j];
			id += 1;
		}
	}
	xadj[graph->N] = id;
}
