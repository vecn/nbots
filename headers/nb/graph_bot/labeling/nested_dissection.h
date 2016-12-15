#ifndef __NB_GRAPH_BOT_LABELING_NESTED_DISSECTION_H__
#define __NB_GRAPH_BOT_LABELING_NESTED_DISSECTION_H__

#include <stdint.h>

void nb_graph_labeling_nd(const nb_graph_t *const graph,
			  uint32_t *perm, uint32_t* iperm);

#endif
