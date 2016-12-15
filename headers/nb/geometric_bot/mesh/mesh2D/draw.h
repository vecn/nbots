#ifndef __NB_GEOMETRIC_BOT_MESH_MESH2D_DRAW_H__
#define __NB_GEOMETRIC_BOT_MESH_MESH2D_DRAW_H__

#include <stdbool.h>

#include "nb/graphics_bot.h"
#include "nb/geometric_bot/mesh/mesh2D.h"
#include "nb/geometric_bot/mesh/mesh2D/info.h"

void nb_mesh2D_export_draw(const nb_mesh2D_t *mesh,
			   const char *filename,
			   int width, int height,
			   nb_mesh2D_entity vals_entity,
			   nb_mesh2D_array_type vals_type,
			   const void *values,
			   bool draw_wires);
void nb_mesh2D_export_level_sets(const nb_mesh2D_t *mesh,
				 const char *filename, int width, int height,
				 const void *field, uint16_t N_ls,
				 bool draw_wires);
void nb_mesh2D_draw_wires(const nb_mesh2D_t *mesh,
			  nb_graphics_context_t *g);
void nb_mesh2D_draw_boundaries(const nb_mesh2D_t *mesh,
			       nb_graphics_context_t *g);
void nb_mesh2D_fill_elems(const nb_mesh2D_t *mesh,
			  nb_graphics_context_t *g);
void nb_mesh2D_fill_elems_field_on_nodes(const nb_mesh2D_t *mesh,
					 nb_graphics_context_t *g,
					 const double *normalized_field,
					 nb_palette_preset palette);
void nb_mesh2D_fill_elems_field_on_elems(const nb_mesh2D_t *mesh,
					 nb_graphics_context_t *g,
					 const double *normalized_field,
					 nb_palette_preset palette);
void nb_mesh2D_fill_elems_classes(const nb_mesh2D_t *mesh,
				  nb_graphics_context_t *g,
				  const uint8_t *class, uint8_t N_colors,
				  const nb_graphics_color_t *colors);
void nb_mesh2D_fill_nodes(const nb_mesh2D_t *mesh,
			  nb_graphics_context_t *g);
void nb_mesh2D_fill_nodes_classes(const nb_mesh2D_t *mesh,
				  nb_graphics_context_t *g,
				  const uint8_t *class, uint8_t N_colors,
				  const nb_graphics_color_t *colors);
void nb_mesh2D_draw_level_set(const nb_mesh2D_t *mesh,
			      nb_graphics_context_t *g,
			      const double *field_on_nodes,
			      double level_set);

#endif
