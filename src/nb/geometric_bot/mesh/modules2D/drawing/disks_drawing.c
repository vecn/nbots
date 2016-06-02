#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>

#include "nb/math_bot.h"
#include "nb/container_bot.h"
#include "nb/geometric_bot.h"
#include "nb/graphics_bot.h"

static void draw_disks(void *draw_ptr, int width, int height,
		       const void *const mshpack_ptr);

void vcn_mshpack_save_png(const vcn_mshpack_t *const mshpack,
			  const char* filename,
			  int width, int height)
{
	nb_drawing_export_png(filename, width, height, draw_disks, mshpack);
}

static void draw_disks(void *draw_ptr, int width, int height,
		       const void *const mshpack_ptr)
{
	const vcn_mshpack_t *const mshpack = mshpack_ptr;

	if (0 == mshpack->N_spheres)
		return;

	/* Compute cam.center and cam.zoom */
	double xmin = mshpack->centers[0] - mshpack->radii[0];
	double xmax = mshpack->centers[0] + mshpack->radii[0];
	double ymin = mshpack->centers[1] - mshpack->radii[0];
	double ymax = mshpack->centers[1] + mshpack->radii[0];
	for (uint32_t i = 1; i < mshpack->N_spheres; i++) {
		if (xmin > mshpack->centers[i * 2] - mshpack->radii[i])
			xmin = mshpack->centers[i * 2] - mshpack->radii[i];
		else if (xmax < mshpack->centers[i * 2] + mshpack->radii[i]) 
			xmax = mshpack->centers[i * 2] + mshpack->radii[i];
		if (ymin > mshpack->centers[i*2+1] - mshpack->radii[i]) 
			ymin = mshpack->centers[i*2+1] - mshpack->radii[i];
		else if (ymax < mshpack->centers[i*2+1] + mshpack->radii[i])
			ymax = mshpack->centers[i*2+1] + mshpack->radii[i];
	}
	camera_t cam;
	cam.center[0] = (xmin+xmax)/2.0;
	cam.center[1] = (ymin+ymax)/2.0;
	cam.zoom = width/(xmax-xmin);
	if (cam.zoom > height/(ymax-ymin))
		cam.zoom = height/(ymax-ymin);
	cam.zoom *= 0.9;


	/* Draw spheres */
	for (uint32_t i = 0; i < mshpack->N_spheres; i++) {
		nb_drawing_set_circle(draw_ptr, &cam,
				      mshpack->centers[i * 2],
				      mshpack->centers[i*2+1],
				      mshpack->radii[i], false);
		
		nb_drawing_set_source_rgba(draw_ptr, 0.1, 0.3, 1.0, 0.5);
		nb_drawing_fill_preserve(draw_ptr);

		nb_drawing_set_line_width(draw_ptr, 0.5);
		nb_drawing_set_source_rgb(draw_ptr, 0, 0, 1);
		nb_drawing_stroke(draw_ptr);
	}
}
