#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>

#include "nb/memory_bot.h"
#include "nb/container_bot.h"
#include "nb/graph_bot/graph.h"

typedef struct {
	uint32_t id;
	double w;
} node_t;

static int8_t node_comparer(const void *n1, const void *n2);
static void extend_adj_1degree(const nb_graph_t *graph,
			       nb_container_t **extended_adj,
			       nb_container_type type,
			       nb_membank_t *membank);
static void set_current_adj(const nb_graph_t *graph,
			    nb_container_t **current_adj,
			    nb_membank_t *membank);
static double get_weight(const nb_graph_t *graph, uint32_t i, uint16_t j);
static void add_1degree_adj(const nb_graph_t *graph,
			    nb_container_t **extended_adj /* NULL */,
			    nb_container_t **current_adj,
			    nb_container_t **new_adj,
			    nb_membank_t *membank);

static void add_nodal_1degree_adj(const nb_graph_t *graph,
				  nb_container_t **extended_adj /* NULL */,
				  nb_container_t **current_adj,
				  nb_container_t **new_adj,
				  uint32_t i, const node_t *node,
				  nb_membank_t *membank);
static bool node_is_not_connected(nb_container_t **extended_adj /* NULL */,
				  nb_container_t **current_adj,
				  nb_container_t **new_adj,
				  uint32_t i, const node_t *node);
static void move_adj(uint32_t N, nb_container_t **dest_adj,
		     nb_container_t **orig_adj);
static void extend_adj(const nb_graph_t *graph,
		       nb_container_t **extended_adj,
		       nb_container_type type,
		       uint8_t N_degrees,
		       nb_membank_t *membank);
static void get_extended_graph(nb_graph_t *graph,
			       nb_container_t **extended_adj,
			       nb_membank_t *membank);
static void realloc_memory(nb_graph_t *graph,
			   nb_container_t **extended_adj);
static uint32_t get_N_total_adj(uint32_t N,
				nb_container_t **extended_adj);


void nb_graph_extend_adj(nb_graph_t *graph, uint8_t N_degrees)
{
	uint32_t N = graph->N;
	nb_container_type type = NB_SORTED;

	uint32_t cnt_size = nb_container_get_memsize(type);
	uint32_t bank_size = nb_membank_get_memsize();
	uint32_t memsize = bank_size +
		N * (sizeof(void*) + cnt_size);
	char *memblock = NB_SOFT_MALLOC(memsize);
	nb_membank_t *membank = (void*) memblock;
	nb_membank_init(membank, sizeof(node_t));
	nb_container_t ** extended_adj = (void*) (memblock + bank_size);
	char *mem = memblock + bank_size + N * sizeof(void*);
	for (uint32_t i = 0; i < N; i++) {
		extended_adj[i] = (void*) mem;
		mem += cnt_size;
		nb_container_init(extended_adj[i], type);
		nb_container_set_comparer(extended_adj[i], node_comparer);
	}
	
	if (1 == N_degrees)
		extend_adj_1degree(graph, extended_adj, type, membank);
	else if (1 < N_degrees)
		extend_adj(graph, extended_adj, type, N_degrees, membank);
	
	get_extended_graph(graph, extended_adj, membank);

	nb_membank_finish(membank);
	for (uint32_t i = 0; i < N; i++)
		nb_container_finish(extended_adj[i]);
	NB_SOFT_FREE(memsize, memblock);
}

static int8_t node_comparer(const void *n1, const void *n2)
{
	const node_t *node1 = n1;
	const node_t *node2 = n2;
	int8_t out;
	if (node1->id < node2->id)
		out = -1;
	else if (node1->id > node2->id)
		out = 1;
	else
		out = 0;
	return out;
}

static void extend_adj_1degree(const nb_graph_t *graph,
			       nb_container_t **extended_adj,
			       nb_container_type type,
			       nb_membank_t *membank)
{
	uint32_t N = graph->N;

	uint32_t cnt_size = nb_container_get_memsize(type);
	uint32_t memsize = N * (sizeof(void*) + cnt_size);
	char *memblock = NB_SOFT_MALLOC(memsize);
	nb_container_t ** new_adj = (void*) memblock;
	char *mem = memblock + N * sizeof(void*);
	for (uint32_t i = 0; i < N; i++) {
		new_adj[i] = (void*) mem;
		mem += cnt_size;
		nb_container_init(new_adj[i], type);
		nb_container_set_comparer(new_adj[i], node_comparer);
	}

	set_current_adj(graph, extended_adj, membank);
	add_1degree_adj(graph, NULL, extended_adj, new_adj, membank);
	move_adj(N, extended_adj, new_adj);

	for (uint32_t i = 0; i < N; i++)
		nb_container_finish(new_adj[i]);
	NB_SOFT_FREE(memsize, memblock);
}

static void set_current_adj(const nb_graph_t *graph,
			    nb_container_t **current_adj,
			    nb_membank_t *membank)
{
	uint32_t N = graph->N;
	for (uint32_t i = 0; i < N; i++) {
		nb_container_t *cnt = current_adj[i];
		for (uint16_t j = 0; j < graph->N_adj[i]; j++) {
			node_t *node = nb_membank_data_calloc(membank);
			node->id = graph->adj[i][j];
			node->w = get_weight(graph, i, j);
			nb_container_insert(cnt, node);
		}
	}
}

static double get_weight(const nb_graph_t *graph, uint32_t i, uint16_t j)
{
	double w;
	if (NULL == graph->wij)
		w = 1.0;
	else
		w = graph->wij[i][j];
	return w;
}

static void add_1degree_adj(const nb_graph_t *graph,
			    nb_container_t **extended_adj /* NULL */,
			    nb_container_t **current_adj,
			    nb_container_t **new_adj,
			    nb_membank_t *membank)
{
	uint32_t N = graph->N;
	nb_iterator_t *iter = alloca(nb_iterator_get_memsize());
	for (uint32_t i = 0; i < N; i++) {
		nb_iterator_init(iter);
		nb_iterator_set_container(iter, current_adj[i]);
		while (nb_iterator_has_more(iter)) {
			const node_t *node = nb_iterator_get_next(iter);
			add_nodal_1degree_adj(graph, extended_adj,
					      current_adj, new_adj,
					      i, node, membank);
		}
		nb_iterator_finish(iter);
	}

}

static void add_nodal_1degree_adj(const nb_graph_t *graph,
				  nb_container_t **extended_adj /* NULL */,
				  nb_container_t **current_adj,
				  nb_container_t **new_adj,
				  uint32_t i, const node_t *node,
				  nb_membank_t *membank)
{
	uint32_t j = node->id;
	double w = node->w;
	for (uint32_t k = 0; k < graph->N_adj[j]; k++) {
		uint32_t id = graph->adj[j][k];
		if (id != i) {
			node_t compare_node;
			compare_node.id = id;
			if (node_is_not_connected(extended_adj, current_adj,
						  new_adj, i, &compare_node)) {
				node_t *new_node =
					nb_membank_data_calloc(membank);
				new_node->id = id;
				new_node->w = w + get_weight(graph, j, k);
				nb_container_insert(new_adj[i], new_node);
			}
		}
	}
}

static bool node_is_not_connected(nb_container_t **extended_adj /* NULL */,
				  nb_container_t **current_adj,
				  nb_container_t **new_adj,
				  uint32_t i, const node_t *node)
{
	bool exist = false;
	if (NULL != extended_adj)
		exist = nb_container_exist(extended_adj[i], node);

	if (!exist)
		exist = nb_container_exist(current_adj[i], node);

	if (!exist)
		exist = nb_container_exist(new_adj[i], node);
	return !exist;
}

static void move_adj(uint32_t N, nb_container_t **dest_adj,
		     nb_container_t **orig_adj)
{
	for (uint32_t i = 0; i < N; i++) {
		while (nb_container_is_not_empty(orig_adj[i])) {
			node_t *node = nb_container_delete_first(orig_adj[i]);
			nb_container_insert(dest_adj[i], node);
		}
	}
}


static void extend_adj(const nb_graph_t *graph,
		       nb_container_t **extended_adj,
		       nb_container_type type,
		       uint8_t N_degrees,
		       nb_membank_t *membank)
{
	uint32_t N = graph->N;

	uint32_t cnt_size = nb_container_get_memsize(type);
	uint32_t memsize = 2 * N * (sizeof(void*) + cnt_size);
	char *memblock = NB_SOFT_MALLOC(memsize);
	nb_container_t ** current_adj = (void*) memblock;
	nb_container_t ** new_adj = (void*) (memblock + N * sizeof(void*));
	char *mem = memblock + 2 * N * sizeof(void*);
	for (uint32_t i = 0; i < N; i++) {
		current_adj[i] = (void*) mem;
		mem += cnt_size;
		nb_container_init(current_adj[i], type);
		nb_container_set_comparer(current_adj[i], node_comparer);

		new_adj[i] = (void*) memblock;
		memblock += cnt_size;
		nb_container_init(new_adj[i], type);
		nb_container_set_comparer(new_adj[i], node_comparer);
	}

	set_current_adj(graph, current_adj, membank);
	for (uint32_t i = 0; i < N_degrees; i++) {
		add_1degree_adj(graph, extended_adj, current_adj,
				new_adj, membank);
		move_adj(N, extended_adj, current_adj);
		move_adj(N, current_adj, new_adj);
	}
	move_adj(N, extended_adj, current_adj);

	for (uint32_t i = 0; i < N; i++) {
		nb_container_finish(current_adj[i]);
		nb_container_finish(new_adj[i]);
	}
	NB_SOFT_FREE(memsize, memblock);
}

static void get_extended_graph(nb_graph_t *graph,
			       nb_container_t **extended_adj,
			       nb_membank_t *membank)
{
	realloc_memory(graph, extended_adj);

	uint32_t N = graph->N;
	for (uint32_t i = 0; i < N; i++) {
		nb_container_t *cnt = extended_adj[i];
		uint32_t id = 0;
		while (nb_container_is_not_empty(cnt)) {
			node_t *node = nb_container_delete_first(cnt);
			graph->adj[i][id] = node->id;
			if (NULL != graph->wij)
				graph->wij[i][id] = node->w;
			id += 1;
			nb_membank_data_free(membank, node);
		}
	}
}

static void realloc_memory(nb_graph_t *graph,
			   nb_container_t **extended_adj)
{
	uint32_t N = graph->N;
	if (NULL != graph->N_adj)
		free(graph->N_adj);
	uint32_t N_total_adj = get_N_total_adj(N, extended_adj);
	uint32_t memsize = N * (sizeof(*(graph->N_adj)) +
				sizeof(*(graph->adj))) +
		N_total_adj * sizeof(**(graph->adj));
	char *memblock = malloc(memsize);
	graph->N_adj = (void*) memblock;
	graph->adj = (void*) (memblock + N * sizeof(*(graph->N_adj)));

	memblock += N * (sizeof(*(graph->N_adj)) + sizeof(*(graph->adj)));
	for (uint32_t i = 0; i < N; i++) {
		graph->adj[i] = (void*) memblock;
		uint32_t N_adj = nb_container_get_length(extended_adj[i]);
		graph->N_adj[i] = N_adj;
		memblock += N_adj * sizeof(**(graph->adj));
	}

	if (NULL != graph->wij) {
		free(graph->wij);
		uint32_t w_memsize = N * sizeof(*(graph->wij)) +
			N_total_adj * sizeof(**(graph->wij));
		char *w_memblock = malloc(w_memsize);
		graph->wij = (void*) w_memblock;
		
		memblock += N * sizeof(*(graph->wij));
		for (uint32_t i = 0; i < N; i++) {
			graph->wij[i] = (void*) memblock;
			uint32_t N_adj = graph->N_adj[i];
			memblock += N_adj * sizeof(**(graph->wij));
		}
	}
}

static uint32_t get_N_total_adj(uint32_t N,
				nb_container_t **extended_adj)
{
	uint32_t N_adj = 0;
	for (uint32_t i = 0; i < N; i++)
		N_adj += nb_container_get_length(extended_adj[i]);
	return N_adj;
}
