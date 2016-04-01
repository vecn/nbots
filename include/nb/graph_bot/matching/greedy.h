#ifndef __NB_GRAPH_BOT_PERFECT_MATCHING_H__
#define __NB_GRAPH_BOT_PERFECT_MATCHING_H__

#include <stdint.h>
#include "nb/graph_bot/graph.h"

void nb_graph_matching_greedy(const nb_graph_t *const graph,
			      int8_t *matched_edges);
