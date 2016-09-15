#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>

#include "nb/memory_bot.h"
#include "nb/container_bot.h"
#include "nb/graph_bot/graph.h"

static void set_current_adj(const nb_graph_t *graph,
			    nb_container_t *current_adj);
static void add_1degree_adj(const nb_graph_t *graph,
			    nb_container_t *current_adj,
			    nb_container_t *new_adj);
static void send_new_to_current_adj(nb_container_t *current_adj,
				    nb_container_t *new_adj);
static void get_extended_graph(nb_graph_t *graph,
			       const nb_container_t *extended_adj);

void nb_graph_extend_adj(nb_graph_t *graph, uint8_t N_degrees)
{
	uint32_t N = graph->N;
	nb_container_type type = NB_SORTED;

	uint32_t cnt_size = nb_container_get_memsize(type);
	uint32_t memsize = N * (sizeof(void*) + cnt_size);
	char *memblock = NB_SOFT_MALLOC(memsize);
	nb_container_t ** extended_adj = (void*) memblock;
	memblock += N * sizeof(void*);
	for (uint32_t i = 0; i < N; i++) {
		extended_adj[i] = (void*) memblock;
		memblock += cnt_size;
		nb_container_init(extended_adj[i], type);
	}

	extend_adj(graph, extended_adj);
	
	get_extended_graph(graph, extended_adj);

	for (uint32_t i = 0; i < N; i++)
		nb_container_finish(extended_adj[i]);
	NB_SOFT_FREE(memsize, memblock);
}

static void extend_adj(const nb_graph_t *graph,
		       nb_container_t **extended_adj)
{
	uint32_t N = graph->N;/* AQUI VOY */
	for (uint32_t i = 0; i < N; i++) {
		set_current_adj(graph, extended_adj);
		for (uint8_t j = 0; j < N_degrees; j++) {
			add_1degree_adj(graph, extended_adj, new_adj);
			send_new_to_current_adj(current_adj, extended_adj);
			send_new_to_current_adj(extended_adj, new_adj);
		}
	}
}

static void set_current_adj(nb_graph_t *graph,
			    const nb_container_t *current_adj)
{
}

static void add_1degree_adj(const nb_graph_t *graph,
			    nb_container_t *current_adj,
			    nb_container_t *extended_adj)
{

}

static void send_new_to_current_adj(nb_container_t *current_adj,
				    nb_container_t *extended_adj);

static void get_extended_graph(nb_graph_t *graph,
			       const nb_container_t **extended_adj)
{

}
