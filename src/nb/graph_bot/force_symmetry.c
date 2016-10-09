#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>

#include "nb/memory_bot.h"
#include "nb/container_bot.h"
#include "nb/graph_bot/graph.h"

static void load_dynamic_sym_graph(nb_container_t **cnt_graph,
				   nb_membank_t *membank,
				   const nb_graph_t *graph);
static void load_sym_graph(uint32_t N,
			   nb_container_t **cnt_graph,
			   nb_membank_t *membank,
			   nb_graph_t *graph);
static void allocate_sym_adj(nb_graph_t *graph,
			     nb_container_t **cnt_graph);
static void set_sym_adj(nb_graph_t *graph,
			nb_container_t **cnt_graph,
			nb_membank_t *membank);

void nb_graph_force_symmetry(nb_graph_t *graph)
{
	nb_container_type cnt_type = NB_SORTED;
	uint32_t cnt_size = nb_container_get_memsize(cnt_type);
	uint32_t N = graph->N;
	uint32_t memsize = N * (sizeof(void*) + cnt_size) +
		nb_membank_get_memsize();
	char *memblock = NB_SOFT_MALLOC(memsize);
	
	nb_container_t **cnt_graph = (void*) memblock;
	for (uint32_t i = 0; i < N; i++) {
		cnt_graph[i] = (void*) (memblock + N * sizeof(void*) +
					i * cnt_size);
		nb_container_init(cnt_graph[i], cnt_type);
	}
	nb_membank_t *membank = (void*) (memblock + N * (sizeof(void*) +
							 cnt_size));
	nb_membank_init(membank, sizeof(uint32_t));

	load_dynamic_sym_graph(cnt_graph, membank, graph);
	nb_graph_clear(graph);
	load_sym_graph(N, cnt_graph, membank, graph);

	for (uint32_t i = 0; i < N; i++)
		nb_container_finish(cnt_graph[i]);

	nb_membank_finish(membank);
	NB_SOFT_FREE(memsize, memblock);
}

static void load_dynamic_sym_graph(nb_container_t **cnt_graph,
				   nb_membank_t *membank,
				   const nb_graph_t *graph)
{
	for (uint32_t i = 0; i < graph->N; i++) {
		uint32_t N_adj = graph->N_adj[i];
		for (uint32_t j = 0; j < N_adj; j++) {
			uint32_t id = graph->adj[i][j];
			uint32_t *adj = nb_membank_allocate_mem(membank);
			*adj = id;
			nb_container_insert(cnt_graph[i], adj);
		}
	}
	for (uint32_t i = 0; i < graph->N; i++) {
		uint32_t N_adj = graph->N_adj[i];
		for (uint32_t j = 0; j < N_adj; j++) {
			uint32_t id = graph->adj[i][j];
			if (NULL == nb_container_exist(cnt_graph[id], &i)) {
				uint32_t *adj =
					nb_membank_allocate_mem(membank);
				*adj = i;
				nb_container_insert(cnt_graph[id], adj);
			}
		}
	}
}

static void load_sym_graph(uint32_t N,
			   nb_container_t **cnt_graph,
			   nb_membank_t *membank,
			   nb_graph_t *graph)
{
	graph->N = N;
	graph->wi = NULL;
	graph->wij = NULL;
	allocate_sym_adj(graph, cnt_graph);
	set_sym_adj(graph, cnt_graph, membank);
}

static void allocate_sym_adj(nb_graph_t *graph,
			     nb_container_t **cnt_graph)
{
	uint32_t total_N_adj = 0;
	for (uint32_t i = 0; i < graph->N; i++)
		total_N_adj += nb_container_get_length(cnt_graph[i]);

	uint32_t memsize_N_adj = graph->N * sizeof(*(graph->N_adj));
	uint32_t memsize_adj = total_N_adj * sizeof(**(graph->adj)) +
		graph->N * sizeof(*(graph->adj));
	char *memblock = nb_allocate_mem(memsize_N_adj + memsize_adj);
	graph->N_adj = (void*) memblock;
	graph->adj = (void*) (memblock + memsize_N_adj);

	memblock += memsize_N_adj + graph->N * sizeof(*(graph->adj));
	for (uint32_t i = 0; i < graph->N; i++) {
		uint32_t N_adj = nb_container_get_length(cnt_graph[i]);
		graph->N_adj[i] = N_adj;
		graph->adj[i] = (void* )memblock;
		memblock += N_adj * sizeof(**(graph->adj));
	}
}

static void set_sym_adj(nb_graph_t *graph,
			nb_container_t **cnt_graph,
			nb_membank_t *membank)
{
	for (uint32_t i = 0; i < graph->N; i++) {
		uint32_t j = 0;
		while (nb_container_is_not_empty(cnt_graph[i])) {
			uint32_t *id = nb_container_delete_first(cnt_graph[i]);
			graph->adj[i][j] = *id;
			j += 1;
			nb_membank_free_mem(membank, id);
		}
	}
}
