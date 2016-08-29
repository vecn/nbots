#ifndef __NB_GEOMETRIC_BOT_MESH_PARTITION_DRAW_H__
#define __NB_GEOMETRIC_BOT_MESH_PARTITION_DRAW_H__

#include <stdbool.h>

#include "nb/geometric_bot/mesh/partition.h"
#include "nb/geometric_bot/mesh/partition/info.h"


void nb_partition_export_draw(const nb_partition_t *part,
			      const char *filename,
			      int width, int height,
			      nb_partition_entity vals_entity,
			      nb_partition_array_type vals_type,
			      const void *values,
			      bool draw_wires);


#endif
