#include <stdlib.h>
#include <stdint.h>

#include "nb/memory_bot.h"
#include "nb/container_bot.h"
#include "nb/geometric_bot/mesh/mesh2D.h"
#include "nb/geometric_bot/mesh/modules2D/graph_generator.h"

#include "../mesh2D_structs.h"

static void allocate_elem_graph(const nb_tessellator2D__t *mesh, nb_graph_t *graph);

void nb_tessellator2D__load_elem_graph(const nb_tessellator2D__t *mesh, nb_graph_t *graph)
{
	uint32_t N_trg = nb_tessellator2D__get_N_trg(mesh);
	if (0 == N_trg) {
		memset(graph, 0, nb_graph_get_memsize());
		goto EXIT;
	}

	graph->N = N_trg;
	graph->wi = NULL;
	graph->wij = NULL;
	allocate_elem_graph(mesh, graph);

	nb_iterator_t* iter = nb_allocate_on_stack(nb_iterator_get_memsize());
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, mesh->ht_trg);
	while (nb_iterator_has_more(iter)) {
		const msh_trg_t* trg = nb_iterator_get_next(iter);
		uint32_t id = trg->id;
		uint32_t cnt = 0;
		if (NULL != trg->t1) {
			uint32_t id2 = trg->t1->id;
			graph->adj[id][cnt] = id2;
			cnt += 1;
		}
		if (NULL != trg->t2) {
			uint32_t id2 = trg->t2->id;
			graph->adj[id][cnt] = id2;
			cnt += 1;
		}
		if (NULL != trg->t3){
			uint32_t id2 = trg->t3->id;
			graph->adj[id][cnt] = id2;
			cnt += 1;
		}
		graph->N_adj[id] = cnt;
	}
	nb_iterator_finish(iter);
EXIT:
	return;
}

static void allocate_elem_graph(const nb_tessellator2D__t *mesh, nb_graph_t *graph)
{
	uint32_t N_trg = nb_tessellator2D__get_N_trg(mesh);
	uint32_t N_total_adj = N_trg * 3;
	uint32_t memsize = N_trg * (sizeof(*(graph->N_adj)) +
				    sizeof(*(graph->adj))) +
		N_total_adj * sizeof(**(graph->adj));

	char *memblock = nb_allocate_mem(memsize);

	graph->N_adj = (void*) memblock;
	graph->adj = (void*) (memblock + N_trg * sizeof(*(graph->N_adj)));

	memblock += N_trg * (sizeof(*(graph->N_adj)) +
				sizeof(*(graph->adj)));
	for (uint32_t i = 0; i < N_trg; i++) {
		graph->adj[i] = (void*) memblock;
		memblock += 3 * sizeof(**(graph->adj));
	}
}
