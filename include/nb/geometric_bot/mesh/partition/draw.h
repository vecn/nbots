#ifndef __NB_GEOMETRIC_BOT_MESH_PARTITION_DRAW_H__
#define __NB_GEOMETRIC_BOT_MESH_PARTITION_DRAW_H__

#include <stdbool.h>

#include "nb/graphics_bot.h"
#include "nb/geometric_bot/mesh/partition.h"
#include "nb/geometric_bot/mesh/partition/info.h"


void nb_partition_export_draw(const nb_partition_t *part,
			      const char *filename,
			      int width, int height,
			      nb_partition_entity vals_entity,
			      nb_partition_array_type vals_type,
			      const void *values,
			      bool draw_wires);
void nb_partition_draw_wires(const nb_partition_t *part,
			     nb_graphics_context_t *g);
void nb_partition_draw_boundaries(const nb_partition_t *part,
				  nb_graphics_context_t *g);
void nb_partition_fill_elems(const nb_partition_t *part,
			     nb_graphics_context_t *g);
void nb_partition_fill_elems_field_on_nodes(const nb_partition_t *part,
					    nb_graphics_context_t *g,
					    const double *normalized_field,
					    nb_palette_preset palette);
void nb_partition_fill_elems_field_on_elems(const nb_partition_t *part,
					    nb_graphics_context_t *g,
					    const double *normalized_field,
					    nb_palette_preset palette);
void nb_partition_fill_elems_classes(const nb_partition_t *part,
				     nb_graphics_context_t *g,
				     const uint8_t *class, uint8_t N_colors,
				     const nb_graphics_color_t *colors);
void nb_partition_fill_nodes(const nb_partition_t *part,
			     nb_graphics_context_t *g);
void nb_partition_fill_nodes_classes(const nb_partition_t *part,
				     nb_graphics_context_t *g,
				     const uint8_t *class, uint8_t N_colors,
				     const nb_graphics_color_t *colors);

#endif
