#include <stdint.h>

#include "nb/memory_bot.h"
#include "nb/graph_bot.h"
#include "nb/eigen_bot/sparse.h"

#include "sparse_struct.h"

void nb_sparse_get_graph(const vcn_sparse_t* A, nb_graph_t *graph)
{
	uint32_t N = A->N;
	graph->N = N;
	uint32_t N_total_adj = vcn_sparse_get_nnz(A) - N;
	uint32_t memsize = N * (sizeof(*(graph->N_adj)) +
				sizeof(*(graph->adj))) +
		N_total_adj * sizeof(**(graph->adj));
	char *memblock = nb_allocate_mem(memsize);

	graph->N_adj = (void*) memblock;
	graph->adj = (void*) (memblock + N * sizeof(*(graph->N_adj)));

	memblock += N * (sizeof(*(graph->N_adj)) + sizeof(*(graph->adj)));
	for (uint32_t i = 0; i < N; i++) {
		uint32_t N_adj = A->rows_size[i] - 1;
		graph->N_adj[i] = N_adj;
		graph->adj[i] = (void*) memblock;
		memblock += N_adj * sizeof(**(graph->adj));
		uint32_t id = 0;
		for (uint32_t j = 0; j < A->rows_size[i]; j++) {
			if (i != A->rows_index[i][j]) {
				graph->adj[i][id] = A->rows_index[i][j];
				id += 1;
			}
		}
	}
}
