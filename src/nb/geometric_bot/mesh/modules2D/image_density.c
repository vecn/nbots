#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#include "nb/math_bot.h"
#include "nb/image_bot.h"
#include "nb/geometric_bot/utils2D.h"
#include "nb/geometric_bot/mesh/mesh2D.h"
#include "nb/geometric_bot/mesh/modules2D/image_density.h"

#include "../mesh2D_structs.h"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

typedef struct {
	const vcn_image_t *img;
	double max_density;
} img_density_data_t;

static void bound_max_density(img_density_data_t *data,
			      double density_volume);

static double img_density(const double x[2],
			  const void *const data_ptr);

void nb_mesh_set_img_density(nb_mesh_t *mesh,
			      const vcn_image_t *const img,
			      double density_volume)
{
	img_density_data_t *data = calloc(1, sizeof(*data));
	data->img = img;
	bound_max_density(data, density_volume);
	nb_mesh_set_density(mesh, img_density, data);
}

static inline void bound_max_density(img_density_data_t *data,
				     double density_volume)
{
	density_volume =  MAX(0.0, density_volume);
	density_volume =  MIN(1.0, density_volume);
	double density_ratio = density_volume * 2.5e-3 +
		(1.0 - density_volume) * 1e-2;
	double max_dim = MAX(vcn_image_get_width(data->img),
			     vcn_image_get_height(data->img));
	double h_min = max_dim * density_ratio;
	data->max_density = 1.0 / h_min;
}

void nb_mesh_clear_img_density(nb_mesh_t *mesh)
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
		return NB_GEOMETRIC_TOL;
	}

	uint8_t luma = vcn_image_get_pixel_luma(data->img, r, c);
	double density = 1.0 - luma / 255.0;
	return data->max_density * density;
}
