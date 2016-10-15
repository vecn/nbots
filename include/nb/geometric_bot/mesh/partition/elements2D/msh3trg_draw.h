#ifndef __NB_GEOMETRIC_BOT_MESH_PARTITION_ELEMENTS2D_MSH3TRG_DRAW_H__
#define __NB_GEOMETRIC_BOT_MESH_PARTITION_ELEMENTS2D_MSH3TRG_DRAW_H__

#include <stdint.h>

#include "nb/graphics_bot.h"

void nb_msh3trg_draw_wires(const void *msh,
			   nb_graphics_context_t *g);
void nb_msh3trg_draw_boundaries(const void *msh,
				nb_graphics_context_t *g);
void nb_msh3trg_fill_elems(const void *msh,
			   nb_graphics_context_t *g);
void nb_msh3trg_fill_elems_field_on_nodes(const void *msh,
					  nb_graphics_context_t *g,
					  const double *normalized_field,
					  nb_palette_preset palette);
void nb_msh3trg_fill_elems_field_on_elems(const void *msh,
					  nb_graphics_context_t *g,
					  const double *normalized_field,
					  nb_palette_preset palette);
void nb_msh3trg_fill_elems_classes(const void *msh,
				   nb_graphics_context_t *g,
				   const uint8_t *class, uint8_t N_colors,
				   const nb_graphics_color_t *colors);
void nb_msh3trg_fill_nodes(const void *msh,
			   nb_graphics_context_t *g);
void nb_msh3trg_fill_nodes_classes(const void *msh,
				   nb_graphics_context_t *g,
				   const uint8_t *class, uint8_t N_colors,
				   const nb_graphics_color_t *colors);

/**
 * @brief Export a PNG image of the partition ans its subdomains.
 * @param[in] part Partition to be displayed.
 * @param[in] k_part Number of partitions.
 * @param[in] part Partition corresponding to the vertex.
 * @param[in] filename Name of the image file.
 * @param[in] width Image width.
 * @param[in] height Image height.
 * @param[in] k_to_draw Partition to draw. Zero to draw all of them.
 * @param[in] scale_partitions Scale of the size of the partitions, in (0,1].
 */
void nb_msh3trg_draw_subdomain(const void *msh3trg_ptr,
			       const char* filename, int width, int height,
			       uint32_t k_part, const uint32_t *const part,
			       uint32_t k_to_draw, double scale_subdomains);

#endif
