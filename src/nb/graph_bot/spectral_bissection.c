#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>

#include "nb/memory_bot.h"
#include "nb/container_bot.h"
#include "nb/solver_bot.h"
#include "nb/graph_bot/graph.h"

static int8_t compare_sb(const void *const A, const void *const B,
			 const void *const data);

static nb_sparse_t* sb_build_laplacian(const nb_graph_t *const graph);
static uint32_t* sb_partition(const nb_graph_t *const graph, uint32_t k,
			  uint32_t N_perm, uint32_t *perm);


static nb_sparse_t* sb_build_laplacian(const nb_graph_t *graph)
{
	nb_sparse_t* B = nb_sparse_create(graph, NULL, 1);
	for (uint32_t i = 0; i < graph->N; i++) {
		nb_sparse_add(B, i, i, 1e-9); /* Damping to the matrix */
		for (uint32_t j = 0; j < graph->N_adj[i]; j++) {
			double wij = 1.0;
			if (graph->wij != NULL)
				wij = graph->wij[i][j];
			nb_sparse_add(B, i, i, wij);
			nb_sparse_add(B, i, graph->adj[i][j], -wij);
		}
	}
	return B;
}

static int8_t compare_sb(const void *const A, const void *const B, 
			 const void *const data)
{
	uint32_t id1 = *((uint32_t*)A);
	uint32_t id2 = *((uint32_t*)B);
	double* fiedler = (double*) data;
	if (fiedler[id1] > fiedler[id2])
		return 1;
	if (fiedler[id1] < fiedler[id2]) 
		return -1;
	return 0;
}

static uint32_t* sb_partition(const nb_graph_t *const __restrict graph, 
			  uint32_t k, uint32_t N_perm, uint32_t *perm)
{
	nb_sparse_t* B = sb_build_laplacian(graph);

	/* Compute eigenvector with second nondecreasing eigenvalue */
	double* eigen_vals = nb_allocate_mem(2 * sizeof(*eigen_vals));
	double** eigen_vecs = nb_allocate_mem(2 * sizeof(*eigen_vecs));
	eigen_vecs[0] = nb_allocate_mem(graph->N * sizeof(*eigen_vecs[0]));
	eigen_vecs[1] = nb_allocate_mem(graph->N * sizeof(*eigen_vecs[1]));

	int iter; /* TEMPORAL: Unused */
	nb_sparse_eigen_ipower(B, NB_SOLVER_LUD,
				2, 0.0, eigen_vecs, eigen_vals,
				&iter, 1e-8, 1);
	nb_sparse_destroy(B);

	uint32_t* p = (uint32_t*) nb_allocate_mem(graph->N * sizeof(uint32_t));
	for (uint32_t i = 0; i < graph->N; i++)
		p[i] = i;

	nb_qsort_wd(p, graph->N, sizeof(*p), compare_sb, eigen_vecs[0]);

	nb_free_mem(eigen_vals);
	nb_free_mem(eigen_vecs[0]);
	nb_free_mem(eigen_vecs[1]);
	nb_free_mem(eigen_vecs);

	/* Recursion */
	uint32_t *k_points1 = NULL;
	uint32_t *k_points2 = NULL;
	if (k > 3) {
		uint32_t mid = N_perm / 2 + N_perm % 2;

		nb_graph_t* subgraph = nb_graph_get_subgraph(graph, mid, p);
		k_points1 = sb_partition(subgraph, k / 2 + k % 2, mid, p);
		nb_free_mem(subgraph);
		
		subgraph = nb_graph_get_subgraph(graph, N_perm / 2, &(p[mid]));
		k_points2 = sb_partition(subgraph, k / 2, N_perm / 2, &(p[mid]));
		nb_free_mem(subgraph);
	}

	/* Set indices in the correct permutation */
	uint32_t* perm_cpy = nb_allocate_mem(N_perm * sizeof(*perm_cpy));
	memcpy(perm_cpy, perm, N_perm * sizeof(*perm_cpy));
	for (uint32_t i = 0; i < graph->N; i++)
		perm[i] = perm_cpy[p[i]];
	nb_free_mem(p);

	uint32_t* k_points = nb_allocate_zero_mem(k * sizeof(*k_points));
	if (k == 2) {
		k_points[0] = N_perm / 2  + N_perm % 2;
		k_points[1] = graph->N ;
	} else if (k == 3) {
		k_points[0] = N_perm / 3 + (N_perm % 3) / 2;
		k_points[1] = 2 * k_points[0];
		k_points[2] = graph->N ;
	} else {
		for (uint32_t i = 0; i < k / 2 + k % 2; i++)
			k_points[i] = k_points1[i];
		for (uint32_t i = 0; i < k / 2; i++) {
			uint32_t mid = N_perm / 2 + N_perm % 2;
			uint32_t id = i + k / 2 + k % 2;
			k_points[id] = mid + k_points2[i];
		}
		nb_free_mem(k_points1);
		nb_free_mem(k_points2);
	}
	return k_points;
}

uint32_t* nb_graph_partition_sb(const nb_graph_t *const graph, uint32_t k)
{
	uint32_t* perm = nb_allocate_mem(graph->N * sizeof(*perm));
	for (uint32_t i = 0; i < graph->N; i++)
		perm[i] = i;

	uint32_t *k_points = sb_partition(graph, k, graph->N, perm);

	uint32_t *part = nb_allocate_mem(graph->N * sizeof(*part));
	uint32_t k_id = 0;
	for (uint32_t i = 0; i < graph->N; i++) {
		if (i >= k_points[k_id])
			k_id += 1;
		part[perm[i]] = k_id;
	}
	nb_free_mem(perm);
	nb_free_mem(k_points);

	return part;
}
