#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#include "vcn/math_bot.h"
#include "vcn/image_bot.h"
#include "vcn/geometric_bot/utils2D.h"
#include "vcn/geometric_bot/mesh/mesh2D.h"
#include "vcn/geometric_bot/mesh/modules2D/image_density.h"

#include "../mesh2D_structs.h"

typedef struct {
	const vcn_image_t *img;
	double max_density;
} img_density_data_t;

static double img_density(const double x[2],
			  const void *const data_ptr);

void vcn_mesh_set_img_density(vcn_mesh_t *mesh,
			      const vcn_image_t *const img,
			      double max_density)
{
	img_density_data_t *data = calloc(1, sizeof(*data));
	data->img = img;
	data->max_density = max_density;
	vcn_mesh_set_density(mesh, img_density, data);
}

void vcn_mesh_clear_img_density(vcn_mesh_t *mesh)
{
	if (NULL != mesh->density_data) {
		free((void*)mesh->density_data);	
		mesh->density_data = NULL;
	}
	mesh->density = NULL;
}

static double img_density(const double x[2],
			  const void *const data_ptr)
{
	const img_density_data_t *const data = data_ptr;

	/* Calculate pixel positions */
	double c_aux = x[0] - vcn_image_get_width(data->img) / 2.0
		+ vcn_image_get_width(data->img) / 2.0;
	double r_aux = x[1] - vcn_image_get_height(data->img) / 2.0
		+ vcn_image_get_height(data->img) / 2.0;
	int c = (int)(c_aux);
	int r = vcn_image_get_height(data->img) - 1 - (int)(r_aux);

	/* Check if the pixel coordinates are outside bounds */
	if (r < 0 || r >= vcn_image_get_height(data->img) ||
	    c < 0 || c >= vcn_image_get_width(data->img)) {
		return VCN_GEOMETRIC_TOL;
	}

	/* Get pixel */
	uint8_t pixel_components[4];
	vcn_image_get_pixel(data->img, r, c, pixel_components);

	double density;
	if (1 == vcn_image_get_N_channels(data->img)) {
		/* Grey */
		density = 1.0 - (1 + pixel_components[0]) / 256.0;
	} else if (2 == vcn_image_get_N_channels(data->img)) {
		/* Grey and alpha */
		double grey = (1 + pixel_components[0]) / 256.0;
		double alpha = (1 + pixel_components[1]) / 256.0;
		density = 1.0 - grey * alpha;
	} else if (3 == vcn_image_get_N_channels(data->img)) {
		/* Red, green and blue */
		double red = (1 + pixel_components[0]) / 256.0;
		double green = (1 + pixel_components[1]) / 256.0;
		double blue = (1 + pixel_components[2]) / 256.0;
		density = 1.0 - (red + green + blue) / 3.0;
	} else {
		/* Red, green, blue and alpha */
		double red = (1 + pixel_components[0]) / 256.0;
		double green = (1 + pixel_components[1]) / 256.0;
		double blue = (1 + pixel_components[2]) / 256.0;
		double alpha = (1 + pixel_components[3]) / 256.0;
		density = 1.0 - alpha * (red + green + blue) / 3.0;
	}
	return data->max_density * density;
}
