#ifndef __VCN_GEOMETRIC_BOT_MESH_EXPORTER_CAIRO_H__
#define __VCN_GEOMETRIC_BOT_MESH_EXPORTER_CAIRO_H__

#include <stdint.h>
#include <stdbool.h>
#include "vcn/geometric_bot.h"

#ifdef __cplusplus
extern "C" {
#endif
	/**
	 * @brief Export a PNG image of the mesh.
	 * @param[in] mesh Mesh to be displayed in the image.
	 * @param[in] filename Name of the PNG file.
	 * @param[in] width Width of the image.
	 * @param[in] height of the image.
	 */
	void vcn_mesh_save_png(const vcn_mesh_t *const mesh,
			       const char* filename,
			       int width, int height);

	void vcn_dewall_save_png(const vcn_mesh_t *const restrict mesh,
				 const char* filename, int width, int height,
				 uint8_t axe, double alpha, uint32_t N,
				 void *vtx_array);
	/**
	 * @brief Export an EPS image of the mesh.
	 * @param[in] mesh Mesh to be displayed in the image.
	 * @param[in] filename Name of the PNG file.
	 * @param[in] width Width of the image.
	 * @param[in] height of the image.
	 */
		void vcn_mesh_save_eps(const vcn_mesh_t *const mesh,
				       const char* filename,
				       int width, int height);

	/**
	 * @brief Export a PNG image of the triangulation.
	 * @param[in] msh3trg Triangulation to be displayed.
	 * @param[in] filename Name of the image file.
	 * @param[in] width Image width.
	 * @param[in] height Image height.
	 * @param[in] rgba_bg Background color, NULL to set white.
	 * @n Color channels: Red, Green, Blue and Alpha.
	 * @param[in] rgba_fg Triangles color, NULL to set default color.
	 * @n Color channels: Red, Green, Blue and Alpha.
	 * @param[in] rgb_edge Edges line color, NULL to set default color.
	 * @n Color channels: Red, Green and Blue.
	 * @param[in] rgb_sgm Segments line color, NULL to set default color.
	 * @n Color channels: Red, Green and Blue.
	 * @param [in] edge_width Line width to draw edges.
	 * @param [in] sgm_width Line width to draw input segments.
	 */
	void vcn_msh3trg_save_png(const vcn_msh3trg_t *const msh3trg,
				  const char* filename, int width, int height,
				  double rgba_bg[4], double rgba_fg[4],
				  double rgb_edge[3], double rgb_sgm[3],
				  double edge_width, double sgm_width);

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
	 * @param[in] rgba_bg Background color, NULL to set white.
	 * @n Color channels: Red, Green, Blue and Alpha.
	 * @param[in] alpha_fg Alpha channel for triangles color.
	 * @n Color channels: Red, Green, Blue and Alpha.
	 * @param[in] rgb_edge Edges line color, NULL to set default color.
	 * @n Color channels: Red, Green and Blue.
	 * @param[in] rgb_sgm Segments line color, NULL to set default color.
	 * @n Color channels: Red, Green and Blue.
	 * @param [in] edge_width Line width to draw edges.
	 * @param [in] sgm_width Line width to draw input segments.
	 */
	void vcn_msh3trg_partition_save_png(const vcn_msh3trg_t *const msh3trg,
					    const char* filename, int width, int height,
					    uint32_t k_part, const uint32_t *const part,
					    uint32_t k_to_draw, double scale_partitions,
					    double rgba_bg[4], double alpha_fg,
					    double rgb_edge[3], double rgb_sgm[3],
					    double edge_width, double sgm_width,
					    bool use_colors);

	/**
	 * @brief Export a PNG image of the Voronoi graph.
	 * @param[in] mshpoly Voronoi graph to be displayed.
	 * @param[in] filename Name of the image file.
	 * @param[in] draw_adjacencies Enable or disable drawing adjacencies.
	 * @param[in] width Image width.
	 * @param[in] height Image height.
	 */
	void vcn_mshpoly_save_png(const vcn_mshpoly_t *const mshpoly, 
				  char* filename,
				  bool draw_adjacencies,
				  int width, int height);

	/**
	 * @brief Export a PNG image of the sphere-pack.
	 * @param[in] mshpack Sphere-pack to be displayed.
	 * @param[in] filename Name of the image file.
	 * @param[in] width Image width.
	 * @param[in] height Image height.
	 */
	void vcn_mshpack_save_png(const vcn_mshpack_t *const mshpack,
				  const char* filename,
				  int width, int height);

#ifdef __cplusplus
}
#endif

#endif
