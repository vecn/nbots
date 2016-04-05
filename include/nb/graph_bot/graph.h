/******************************************************************************
 *   Graph's Bot: Fast graph theory tools.                                    *
 *   2011-2015 Victor Eduardo Cardoso Nungaray                                *
 *   Twitter: @victore_cardoso                                                *
 *   email: victorc@cimat.mx                                                  *
 ******************************************************************************/

/**
 * @file graph_bot.h
 * @brief The graph bot is a set of numerical procedures to handle graphs,
 * exahistively used in meshing procedures, parallel computations and machine
 * learning techniques.
 * @author Victor Eduardo Cardoso Nungaray
 * @n victorc@@cimat.mx
 * @n <a href="https://twitter.com/victore_cardoso"> @@victore_cardoso </a>
 * @date 10 August 2015
 *
 * @mainpage Graph's Bot
 * A graph operations for numerical analysis.
 */

#ifndef __NB_GRAPH_BOT_GRAPH_H__
#define __NB_GRAPH_BOT_GRAPH_H__

#include <stdint.h>

typedef struct {
  uint32_t N;
  uint32_t *N_adj;
  uint32_t **adj;
  double *wi;   /* NULL for equal weights on vertices */
  double **wij; /* NULL for equal weights on edges */
} nb_graph_t;

typedef nb_graph_t vcn_graph_t; /* Deprecated */

vcn_graph_t* vcn_graph_create(void);

/**
 * @brief Destroy graph.
 * @param[in] graph Graph to be destroyed.
 */
void vcn_graph_destroy(vcn_graph_t *graph);

void nb_graph_init_vtx_weights(nb_graph_t *graph);
void nb_graph_init_edge_weights(nb_graph_t *graph);
void nb_graph_finish_vtx_weights(nb_graph_t *graph);
void nb_graph_finish_edge_weights(nb_graph_t *graph);

/**
 * @brief Get a subgraph from the graph.
 * @param[in] graph Super graph containing the output subgraph.
 * @param[in] N_nodes Number of nodes in the subgraph.
 * @param[in] nodes ID of nodes of the subgraph.
 * @return Subgraph.
 */
vcn_graph_t* vcn_graph_get_subgraph(const vcn_graph_t *const graph,
				    uint32_t N_nodes, uint32_t *nodes);

/**
 * @brief Labeling used to minimize entries in LU decomposition
 * based on the AMD algorithm.
 * @param[in] graph Input graph to be labeled.
 * @param[out] iperm Inverse permutation, array mapping output IDs
 * with input IDs. NULL if not required.
 * @return The permutation, an array mapping input IDs with output IDs.
 *
 * @see P. R. Amestoy, T. A. Davis and I. S. Duff. <b>An Approximate
 * Minimum Degree Ordering Algorithm.</b> Journal of Matrix Analysis
 * and Applications (Vol. 17 1996), SIAM, pages 886-905.
 */
uint32_t* vcn_graph_labeling_amd(const vcn_graph_t *const graph, uint32_t* iperm);

/**
 * @brief Multiple Minimum Degree based on Quotient Graphs. The output
 * labeling reduces the fill-in on LU decompositions.
 * @param[in] graph Input graph to be labeled.
 * @param[out] iperm Inverse permutation, array mapping output IDs
 * with input IDs. NULL if not required.
 * @return The permutation, an array mapping input IDs with output IDs.
 *
 * @see A. George and J. W.H. Liu. <b>The evolution of the minimum degree
 * ordering algorithm.</b> SIAM Review (Vol. 31 1989), pages 1-19.
 */
uint32_t* vcn_graph_labeling_mmd(const vcn_graph_t *const graph, uint32_t* iperm);

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
uint32_t* vcn_graph_partition_sb(const vcn_graph_t *const graph, uint32_t k);

#endif
