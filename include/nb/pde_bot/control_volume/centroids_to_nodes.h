#ifndef __NB_PDE_BOT_CONTROL_VOLUME_CENTROIDS_TO_NODES_H__
#define __NB_PDE_BOT_CONTROL_VOLUME_CENTROIDS_TO_NODES_H__

#include <stdint.h>

#include "nb/geometric_bot.h"

int nb_cvfa_interpolate_from_centroids_to_nodes
			(const nb_partition_t *const part,
			 uint32_t N_comp,
			 const double* cen_values,
			 double* nodal_values /* Output */);

#endif
