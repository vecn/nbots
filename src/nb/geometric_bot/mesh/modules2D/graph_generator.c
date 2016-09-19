#include <stdlib.h>
#include <stdint.h>

#include "nb/container_bot.h"
#include "nb/geometric_bot/mesh/mesh2D.h"
#include "nb/geometric_bot/mesh/modules2D/graph_generator.h"

#include "../mesh2D_structs.h"

nb_graph_t* vcn_mesh_create_vtx_graph(const vcn_mesh_t *const restrict mesh)
{
	nb_graph_t* graph = calloc(1, sizeof(*graph));

	graph->N = vcn_bins2D_get_length(mesh->ug_vtx);
	graph->N_adj = calloc(graph->N, sizeof(*(graph->N_adj)));
	graph->adj = malloc(graph->N * sizeof(*(graph->adj)));
	nb_iterator_t* iter = nb_iterator_create();
	nb_iterator_set_container(iter, mesh->ht_edge);
	while (nb_iterator_has_more(iter)) {
		msh_edge_t* edge = (msh_edge_t*)nb_iterator_get_next(iter);
		uint32_t idx1 = mvtx_get_id(edge->v1);
		uint32_t idx2 = mvtx_get_id(edge->v2);
		graph->N_adj[idx1] += 1;
		graph->N_adj[idx2] += 1;
	}

	for (uint32_t i = 0; i < graph->N; i++)
		graph->adj[i] = malloc(graph->N_adj[i] * sizeof(*(graph->adj[i])));

	uint32_t* adj_next_idx = calloc(graph->N, sizeof(*adj_next_idx));
	nb_iterator_restart(iter);
	while (nb_iterator_has_more(iter)) {
		msh_edge_t* edge = (msh_edge_t*)nb_iterator_get_next(iter);
		uint32_t idx1 = mvtx_get_id(edge->v1);
		uint32_t idx2 = mvtx_get_id(edge->v2);
		graph->adj[idx1][adj_next_idx[idx1]] = idx2;
		graph->adj[idx2][adj_next_idx[idx2]] = idx1;
		adj_next_idx[idx1] += 1;
		adj_next_idx[idx2] += 1;
	}
	nb_iterator_destroy(iter);
	free(adj_next_idx);

	return graph;
}

nb_graph_t* vcn_mesh_create_elem_graph(const vcn_mesh_t *const restrict mesh)
{
	nb_graph_t* graph = calloc(1, sizeof(*graph));
	graph->N = nb_container_get_length(mesh->ht_trg);
	graph->N_adj = calloc(graph->N, sizeof(*(graph->N_adj)));
	graph->adj = malloc(graph->N * sizeof(*(graph->adj)));
	nb_iterator_t* iter = nb_iterator_create();
	nb_iterator_set_container(iter, mesh->ht_trg);
	while (nb_iterator_has_more(iter)) {
		msh_trg_t* trg = (msh_trg_t*)nb_iterator_get_next(iter);
		uint32_t id = trg->id;
		if (NULL != trg->t1)
			graph->N_adj[id] += 1;
		if (NULL != trg->t2)
			graph->N_adj[id] += 1;
		if (NULL != trg->t3)
			graph->N_adj[id] += 1;
	}

	nb_iterator_restart(iter);
	while (nb_iterator_has_more(iter)) {
		msh_trg_t* trg = (msh_trg_t*)nb_iterator_get_next(iter);
		uint32_t id = trg->id;
		graph->adj[id] = malloc(graph->N_adj[id] * sizeof(*(graph->adj[id])));
		int cnt = 0;
		if (NULL != trg->t1) {
			uint32_t id2 = trg->t1->id;
			graph->adj[id][cnt++] = id2;
		}
		if (NULL != trg->t2) {
			uint32_t id2 = trg->t2->id;
			graph->adj[id][cnt++] = id2;
		}
		if (NULL != trg->t3){
			uint32_t id2 = trg->t3->id;
			graph->adj[id][cnt++] = id2;
		}
	}
	nb_iterator_destroy(iter);
	return graph;
}
