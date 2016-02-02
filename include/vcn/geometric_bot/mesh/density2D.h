#ifndef __VCN_DENSITY_2D_H__
#define __VCN_DENSITY_2D_H__

#include <stdint.h>

/**
 * @brief Set as density function in vcn_mesh_refine() to create the
 * Constrained Delaunay Triangulation of the domain, without inserting
 * new vertices.
 */
#define VCN_DENSITY_CDT ((void*)1)

/**
 * @brief Set as density function in vcn_mesh_refine() to create a mesh by
 * defining the max size of edges in the output.
 */
#define VCN_DENSITY_MAX ((void*)2)

/**
 * @brief Set as density function in vcn_mesh_refine() to utilize an image
 * as density function.
 */
#define VCN_DENSITY_IMG ((void*)3)

#ifdef __cplusplus
extern "C" {
#endif

	/**
	 * @brief Built-in image used to calculate the density of a 
	 * given point.
	 */
	typedef struct vcn_density_img_s vcn_density_img_t;
  

	/**
	 * @brief Read the image to be used in the density calculation.
	 * @param[in] filename File name of the image. Formats supported:
	 * PNG, JPEG, BMP, PSD (composite view only, no extra channels), GIF,
	 * HDR, PIC and PNM (PPM and PGM binary only). 
	 * @param[in] scale Scale of the image. To keep the same size, set a scale
	 * of 1.
	 * @param[in[ xdisp Horizontal displacement of the image. This parameter
	 * could be used for centering the image along the X axis.
	 * @param[in] ydisp Vertical displacement of the image. This parameter could
	 * be used for centering the image along the Y axis.
	 * @param[in] max_density The highest density representable in the image. The
	 * max_density corresponds to black and the lowest density, zero, corresponds
	 * to white.
	 * @return The image data required when the density parameter is set to 2
	 * (density given by an image).
	 */
	vcn_density_img_t* vcn_density_img_create(const char* filename,
						  double scale,
						  double xdisp,
						  double ydisp,
						  double max_density);

	/**
	 * @brief Returns the image's width.
	 * @param[in] data image's data.
	 */
	uint32_t vcn_density_img_get_width(const vcn_density_img_t *const data);

	/**
	 * @brief Returns the image's height.
	 * @param[in] data image's data.
	 */
	uint32_t vcn_density_img_get_height
				(const vcn_density_img_t *const data);
	double vcn_density_img_get_density(const double x[2],
					   const void *const density2D_ptr);
	/**
	 * @brief Destroys the density image data.
	 * @param[in] data image's data to be destroyed (memory free).
	 */
	void vcn_density_img_destroy(vcn_density_img_t* data);

#ifdef __cplusplus
}
#endif

#endif
