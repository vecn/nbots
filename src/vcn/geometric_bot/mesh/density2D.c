/******************************************************************************
 *   Geometric Bot: Geometric tesselations for numerical analysis.            *
 *   2011-2015 Victor Eduardo Cardoso Nungaray                                *
 *   Twitter: @victore_cardoso                                                *
 *   email: victorc@cimat.mx                                                  *
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#include "vcn/math_bot.h"
#include "vcn/geometric_bot/utils2D.h"
#include "vcn/geometric_bot/mesh/density2D.h"

#define _VCN_INPUT_VTX ((void*)0x1)
#define _VCN_SUBSEGMENT_VTX ((void*)0x2)
#define _VCN_CC_SHELL_UNIT (1e-3)
#define _VCN_MAX_GRADING_RATIO (27.0)
#define _VCN_MAX_LH_TOLERATED (1.5)

#define STB_IMAGE_IMPLEMENTATION
#include "ext/stb_image.h"
   
/******************** Private structures *****************************/

struct vcn_density_img_s {
	uint32_t width;
	uint32_t height;
	uint32_t comp_x_pixel;
	double scale;
	double xdisp, ydisp;
	double max_density;
	uint8_t* pixels;
};

vcn_density_img_t* vcn_density_img_create(const char* filename,
					  double scale,
					  double xdisp,
					  double ydisp,
					  double max_density)
{
	/* Read image */
	int w, h, n;
	uint8_t *const restrict pixels = stbi_load(filename, &w, &h, &n, 0);

	if (pixels == NULL) 
		return NULL;

	/* Allocate structure */
	vcn_density_img_t* data = malloc(sizeof(vcn_density_img_t));
	data->pixels = pixels;
	data->width = w;
	data->height = h;
	data->comp_x_pixel = n;
	data->scale = scale;
	data->xdisp = xdisp;
	data->ydisp = ydisp;
	data->max_density = max_density;

	/* Return structure */
	return data;
}

inline uint32_t vcn_density_img_get_width(const vcn_density_img_t *const data)
{
	return data->width;
}

inline uint32_t vcn_density_img_get_height(const vcn_density_img_t *const data)
{
	return data->height;
}

double vcn_density_img_get_density(const void *const density2D_ptr,
				   const double x[2])
{
	const vcn_density_img_t *const density2D = density2D_ptr;

	/* Calculate pixel positions (Centered-Scaling and displacement) */
	double c_aux = (x[0] - density2D->width / 2.0)  /
		density2D->scale + density2D->width / 2.0;
	double r_aux = (x[1] - density2D->height / 2.0) / 
		density2D->scale + density2D->height / 2.0;
	int c = (int)(c_aux - density2D->xdisp);
	int r = density2D->height - 1 - (int)(r_aux - density2D->ydisp);

	/* Check if the pixel coordinates are outside bounds */
	if (r < 0 || r >= density2D->height || c < 0 || c >= density2D->width)
		return VCN_GEOMETRIC_TOL;

	/* Get pixel */
	uint8_t pixel_components[4];
	for (uint8_t i = 0; i < density2D->comp_x_pixel; i++) {
		uint32_t w = density2D->width * density2D->comp_x_pixel;
		uint32_t col = c * density2D->comp_x_pixel;
		pixel_components[i] = density2D->pixels[r * w + col + i];
	}

	double density;
	if (density2D->comp_x_pixel == 1) {
		/* Grey */
		density = density2D->max_density * 
			(1.0 - (1 + pixel_components[0]) / 256.0);
	} else if (density2D->comp_x_pixel == 2) {
		/* Grey and alpha */
		double grey = (1 + pixel_components[0]) / 256.0;
		double alpha = (1 + pixel_components[1]) / 256.0;
		density = density2D->max_density * (1.0 - grey * alpha);
	} else if(density2D->comp_x_pixel == 3) {
		/* Red, green and blue */
		double red = (1 + pixel_components[0]) / 256.0;
		double green = (1 + pixel_components[1]) / 256.0;
		double blue = (1 + pixel_components[2]) / 256.0;
		density = density2D->max_density * (1.0 - (red + green + blue) / 3.0);
	} else {
		/* Red, green, blue and alpha */
		double red = (1 + pixel_components[0]) / 256.0;
		double green = (1 + pixel_components[1]) / 256.0;
		double blue = (1 + pixel_components[2]) / 256.0;
		double alpha = (1 + pixel_components[3]) / 256.0;
		density =  density2D->max_density *
			(1.0 - alpha * (red + green + blue) / 3.0);
	}
	return density;
}

void vcn_density_img_destroy(vcn_density_img_t* data)
{
	stbi_image_free(data->pixels);
	free(data);
}
