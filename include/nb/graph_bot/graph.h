/**
 * @file graph_bot.h
 * @brief The graph bot is a set of numerical procedures to handle graphs,
 * exahistively used in meshing procedures, parallel computations and machine
 * learning techniques.
 * @author Victor Eduardo Cardoso Nungaray
 * @n victorc@@cimat.mx
 * @n <a href="https://twitter.com/victore_cardoso"> @@victore_cardoso </a>
 *
 * @mainpage Graph's Bot
 * A graph operations for numerical analysis.
 */

#ifndef __NB_GRAPH_BOT_GRAPH_H__
#define __NB_GRAPH_BOT_GRAPH_H__

#include <stdint.h>

typedef enum {
	NB_LABELING_MMD,
	NB_LABELING_AMD,
	NB_LABELING_ND
} nb_labeling_algorithm;

typedef struct {
	uint32_t N;
	uint32_t *N_adj;
	uint32_t **adj;
	double *wi;   /* NULL for equal weights on vertices */
	double **wij; /* NULL for equal weights on edges */
} nb_graph_t;

uint32_t nb_graph_get_memsize(void);
void nb_graph_init(nb_graph_t *graph);
void nb_graph_copy(void *graph, const void *graph_src);
void nb_graph_finish(nb_graph_t *graph);
void nb_graph_clear(nb_graph_t *graph);

void nb_graph_init_vtx_weights(nb_graph_t *graph);
void nb_graph_init_edge_weights(nb_graph_t *graph);
void nb_graph_finish_vtx_weights(nb_graph_t *graph);
void nb_graph_finish_edge_weights(nb_graph_t *graph);

void nb_graph_force_symmetry(nb_graph_t *graph);
void nb_graph_extend_adj(nb_graph_t *graph, uint8_t N_degrees);

uint32_t nb_graph_find_N_intersected_adj(const nb_graph_t *graph,
					uint32_t id1, uint32_t id2);

/**
 * @brief Get a subgraph from the graph.
 * @param[in] graph Super graph containing the output subgraph.
 * @param[in] N_nodes Number of nodes in the subgraph.
 * @param[in] nodes ID of nodes of the subgraph.
 * @return Subgraph.
 */
nb_graph_t* nb_graph_get_subgraph(const nb_graph_t *const graph,
				    uint32_t N_nodes, uint32_t *nodes);

void nb_graph_labeling(const nb_graph_t *const graph,
		       uint32_t *perm, uint32_t* iperm,
		       nb_labeling_algorithm alg);

/**
 * @brief Graph partition using Spectral bisection.
 * @param[in] graph Graph to be partitioned.
 * @return Arrays of IDs corresponding to the partition of each node, starting
 * from zero.
 *
 * @see B. Hendrickson and R. Leland. <b>An improved spectral graph partitioning 
 * algorithm for mapping parallel computations.</b> SIAM J. Sci. Comput.
 * (Vol 16 1995), pages 452-469.
 */
uint32_t* nb_graph_partition_sb(const nb_graph_t *const graph, uint32_t k);

#endif
