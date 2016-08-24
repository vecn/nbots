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
		       const void * mshpack);

void nb_mshpack_draw(const void *mshpack, const char* filename,
		     int width, int height)
{
	nb_graphics_export(filename, width, height, draw_disks, mshpack);
}

static void draw_disks(nb_graphics_context_t *g, int width, int height,
		       const void *mshpack)
{
	if (0 < nb_mshpack_get_N_elems(mshpack)) {
		double box[4];
		nb_mshpack_get_enveloping_box(mshpack, box);

		nb_graphics_enable_camera(g);
		nb_graphics_camera_t* cam = nb_graphics_get_camera(g);
		nb_graphics_cam_fit_box(cam, box, width, height);

		uint32_t N_elems = nb_mshpack_get_N_elems(mshpack);
		for (uint32_t i = 0; i < N_elems; i++) {
			double x = nb_mshpack_get_x_elem(mshpack, i);
			double y = nb_mshpack_get_y_elem(mshpack, i);
			double r = nb_mshpack_elem_get_radii(mshpack, i);
			nb_graphics_set_circle(g, x, y, r);
		
			nb_graphics_set_source_rgba(g, 25, 75, 255, 128);
			nb_graphics_fill_preserve(g);

			nb_graphics_set_line_width(g, 0.5);
			nb_graphics_set_source(g, NB_BLUE);
			nb_graphics_stroke(g);
		}
	}
}
