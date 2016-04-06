#ifndef __NB_GRAPH_BOT_MATCHING_GREEDY_H__
#define __NB_GRAPH_BOT_MATCHING_GREEDY_H__

#include <stdint.h>
#include "nb/graph_bot/graph.h"

void nb_graph_matching_greedy(const nb_graph_t *const graph,
			      uint32_t *nodal_match);
