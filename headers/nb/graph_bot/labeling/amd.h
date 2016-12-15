#ifndef __NB_GRAPH_BOT_LABELING_AMD_H__
#define __NB_GRAPH_BOT_LABELING_AMD_H__

#include <stdint.h>

void nb_graph_labeling_amd(const nb_graph_t *const graph,
			   uint32_t *perm, uint32_t* iperm);

void nb_graph_labeling_mmd(const nb_graph_t *const graph,
			   uint32_t *perm, uint32_t* iperm);

#endif
