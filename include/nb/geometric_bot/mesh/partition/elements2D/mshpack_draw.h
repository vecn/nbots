#ifndef __NB_GEOMETRIC_BOT_MESH_PARTITION_ELEMENTS2D_MSHPACK_DRAW_H__
#define __NB_GEOMETRIC_BOT_MESH_PARTITION_ELEMENTS2D_MSHPACK_DRAW_H__

#include <stdint.h>

#include "nb/graphics_bot.h"

void nb_mshpack_draw_wires(const void *msh,
			   nb_graphics_context_t *g);
void nb_mshpack_draw_boundaries(const void *msh,
				nb_graphics_context_t *g);
void nb_mshpack_fill_elems(const void *msh,
			   nb_graphics_context_t *g);
void nb_mshpack_fill_elems_field_on_nodes(const void *msh,
					  nb_graphics_context_t *g,
					  const double *normalized_field,
					  nb_palette_preset palette);
void nb_mshpack_fill_elems_field_on_elems(const void *msh,
					  nb_graphics_context_t *g,
					  const double *normalized_field,
					  nb_palette_preset palette);
void nb_mshpack_fill_elems_classes(const void *msh,
				   nb_graphics_context_t *g,
				   const uint8_t *class, uint8_t N_colors,
				   const nb_graphics_color_t *colors);
void nb_mshpack_fill_nodes(const void *msh,
			   nb_graphics_context_t *g);
void nb_mshpack_fill_nodes_classes(const void *msh,
				   nb_graphics_context_t *g,
				   const uint8_t *class, uint8_t N_colors,
				   const nb_graphics_color_t *colors);

#endif
