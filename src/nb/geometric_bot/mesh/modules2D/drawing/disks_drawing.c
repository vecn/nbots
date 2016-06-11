#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>

#include "nb/math_bot.h"
#include "nb/container_bot.h"
#include "nb/geometric_bot.h"
#include "nb/graphics_bot.h"

static void draw_disks(nb_graphics_context_t *g, int width, int height,
		       const void *const mshpack_ptr);

void vcn_mshpack_draw(const vcn_mshpack_t *const mshpack,
			  const char* filename,
			  int width, int height)
{
	nb_graphics_export(filename, width, height, draw_disks, mshpack);
}

static void draw_disks(nb_graphics_context_t *g, int width, int height,
		       const void *const mshpack_ptr)
{
	const vcn_mshpack_t *const mshpack = mshpack_ptr;

	if (0 == mshpack->N_spheres)
		return;

	/* Compute cam.center and cam.zoom */
	double box[4];
	box[0] = mshpack->centers[0] - mshpack->radii[0];
	box[1] = mshpack->centers[1] - mshpack->radii[0];
	box[2] = mshpack->centers[0] + mshpack->radii[0];
	box[3] = mshpack->centers[1] + mshpack->radii[0];
	for (uint32_t i = 1; i < mshpack->N_spheres; i++) {
		if (box[0] > mshpack->centers[i * 2] - mshpack->radii[i])
			box[0] = mshpack->centers[i * 2] - mshpack->radii[i];
		else if (box[2] < mshpack->centers[i * 2] + mshpack->radii[i]) 
			box[2] = mshpack->centers[i * 2] + mshpack->radii[i];
		if (box[1] > mshpack->centers[i*2+1] - mshpack->radii[i]) 
			box[1] = mshpack->centers[i*2+1] - mshpack->radii[i];
		else if (box[3] < mshpack->centers[i*2+1] + mshpack->radii[i])
			box[3] = mshpack->centers[i*2+1] + mshpack->radii[i];
	}

	nb_graphics_enable_camera(g);
	camera_t* cam = nb_graphics_get_camera(g);
	nb_graphics_cam_fit_box(cam, box, width, height);

	/* Draw spheres */
	for (uint32_t i = 0; i < mshpack->N_spheres; i++) {
		nb_graphics_set_circle(g,
				       mshpack->centers[i * 2],
				       mshpack->centers[i*2+1],
				       mshpack->radii[i]);
		
		nb_graphics_set_source_rgba(g, 25, 75, 255, 128);
		nb_graphics_fill_preserve(g);

		nb_graphics_set_line_width(g, 0.5);
		nb_graphics_set_source(g, NB_BLUE);
		nb_graphics_stroke(g);
	}
}
