#ifndef __NB_GEOMETRIC_BOT_MESH_MODULES2D_DRAWING_H__
#define __NB_GEOMETRIC_BOT_MESH_MODULES2D_DRAWING_H__

#include <stdint.h>
#include <stdbool.h>
#include "nb/geometric_bot/mesh/elements2D/triangles_struct.h"
#include "nb/geometric_bot/mesh/elements2D/polygons_struct.h"
#include "nb/geometric_bot/mesh/elements2D/quad_struct.h"
#include "nb/geometric_bot/mesh/elements2D/disks_struct.h"

void vcn_mesh_draw(const vcn_mesh_t *const mesh,
		   const char* filename,
		   int width, int height);

void vcn_dewall_draw(const vcn_mesh_t *const restrict mesh,
		     const char* filename, int width, int height,
		     uint8_t axe, double alpha, uint32_t N,
		     void *vtx_array);

void vcn_msh3trg_draw(const vcn_msh3trg_t *const msh3trg,
		      const char* filename, int width, int height);

/**
 * @brief Export a PNG image of the triangulation and its partitions.
 * @param[in] msh3trg Triangulation to be displayed.
 * @param[in] k_part Number of partitions.
 * @param[in] part Partition corresponding to the vertex.
 * @param[in] filename Name of the image file.
 * @param[in] width Image width.
 * @param[in] height Image height.
 * @param[in] k_to_draw Partition to draw. Zero to draw all of them.
 * @param[in] scale_partitions Scale of the size of the partitions, in (0,1].
 */
void vcn_msh3trg_partition_draw(const vcn_msh3trg_t *const msh3trg,
				const char* filename, int width, int height,
				uint32_t k_part, const uint32_t *const part,
				uint32_t k_to_draw, double scale_partitions);

void nb_mshquad_draw(const nb_mshquad_t *const quad,
		     const char* filename, int width, int height);

void nb_mshpoly_draw(const nb_mshpoly_t *const poly,
			   const char* filename, int width, int height);

void vcn_mshpack_draw(const vcn_mshpack_t *const mshpack,
		      const char* filename,
		      int width, int height);

#endif
