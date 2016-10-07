#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#ifndef __NB_PDE_BOT_CONTROL_VOLUME_INTEGRATION_MESH_H__
#define __NB_PDE_BOT_CONTROL_VOLUME_INTEGRATION_MESH_H__

#include <stdint.h>

#include "nb/geometric_bot.h"

uint32_t nb_cvfa_get_integration_mesh_memsize(void);
void nb_cvfa_init_integration_mesh(nb_partition_t *intmsh);
void nb_cvfa_load_integration_mesh(const nb_partition_t *part,
				   nb_partition_t *intmsh);

void nb_cvfa_correlate_partition_and_integration_mesh
					(const nb_partition_t *part,
					 const nb_partition_t *intmsh,
					 nb_graph_t *trg_x_vol);

void nb_cvfa_get_adj_graph(const nb_partition_t *intmsh,
			   const nb_graph_t *trg_x_vol,
			   nb_graph_t *graph);

#endif
