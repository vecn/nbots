#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>

#include "nb/math_bot.h"
#include "nb/container_bot.h"
#include "nb/geometric_bot.h"
#include "nb/graphics_bot.h"

static void draw_mesh(void *draw_ptr, int width, int height,
		      const void *const quad_ptr);
static void draw_quads(void *draw_ptr, const nb_mshquad_t *const quad,
		       const camera_t *cam);
static void draw_input_sgm(void *draw_ptr, const nb_mshquad_t *const quad,
			   const camera_t *cam);
static void draw_input_vtx(void *draw_ptr, const nb_mshquad_t *const quad,
			   const camera_t *cam);

void nb_mshquad_export_png(const nb_mshquad_t *const quad,
			   const char* filename, int width, int height)
{
	nb_drawing_export_png(filename, width, height, draw_mesh, quad);
}

void nb_mshquad_export_eps(const nb_mshquad_t *const quad,
			   const char* filename, int width, int height)
{
	nb_drawing_export_eps(filename, width, height, draw_mesh, quad);
}

static void draw_mesh(void *draw_ptr, int width, int height,
		      const void *const quad_ptr)
{
	const nb_mshquad_t *const quad = quad_ptr;

	double box[4];
	vcn_utils2D_get_enveloping_box_from_subset(quad->N_vtx, quad->vtx,
						   quad->nod,
						   2 * sizeof(*(quad->nod)),
						   vcn_utils2D_get_x_from_darray,
						   vcn_utils2D_get_y_from_darray,
						   box);
	camera_t cam;
	nb_drawing_utils_set_center_and_zoom(&cam, box, width, height);
	draw_quads(draw_ptr, quad, &cam);
	draw_input_sgm(draw_ptr, quad, &cam);
	draw_input_vtx(draw_ptr, quad, &cam);
}

static void draw_quads(void *draw_ptr, const nb_mshquad_t *const quad,
		       const camera_t *cam)
{
	nb_drawing_set_line_width(draw_ptr, 0.5);
	for (uint32_t i = 0; i < quad->N_elems; i++) {
		uint32_t n1 = quad->adj[i * 4];
		uint32_t n2 = quad->adj[i*4+1];
		uint32_t n3 = quad->adj[i*4+2];
		uint32_t n4 = quad->adj[i*4+3];
		nb_drawing_move_to(draw_ptr, cam,
				   quad->nod[n1 * 2],
				   quad->nod[n1*2+1]);
		nb_drawing_line_to(draw_ptr, cam,
				   quad->nod[n2 * 2],
				   quad->nod[n2*2+1]);
		nb_drawing_line_to(draw_ptr, cam,
				   quad->nod[n3 * 2],
				   quad->nod[n3*2+1]);
		if (n4 < quad->N_nod)
			nb_drawing_line_to(draw_ptr, cam,
					   quad->nod[n4 * 2],
					   quad->nod[n4*2+1]);
		nb_drawing_close_path(draw_ptr);

		nb_drawing_set_source_rgba(draw_ptr, 0.1, 0.3, 1.0, 0.5);
		nb_drawing_fill_preserve(draw_ptr);

		nb_drawing_set_source_rgb(draw_ptr, 0.0, 0.0, 1.0);
		nb_drawing_stroke(draw_ptr);
	}
}

static void draw_input_sgm(void *draw_ptr, const nb_mshquad_t *const quad,
			   const camera_t *cam)
{
	nb_drawing_set_line_width(draw_ptr, 1.0);
	nb_drawing_set_source_rgb(draw_ptr, 0.9, 0.2, 0.4);
	for (uint32_t i = 0; i < quad->N_sgm; i++) {
		if (0 < quad->N_nod_x_sgm[i]) {
			uint32_t nj = quad->nod_x_sgm[i][0];
			nb_drawing_move_to(draw_ptr, cam,
					   quad->nod[nj * 2],
					   quad->nod[nj*2+1]);
			for (uint32_t j = 1; j < quad->N_nod_x_sgm[i]; j++) {
				nj = quad->nod_x_sgm[i][j];
				nb_drawing_line_to(draw_ptr, cam,
						   quad->nod[nj * 2],
						   quad->nod[nj*2+1]);
			}
			nb_drawing_stroke(draw_ptr);
		}
	}
}

static void draw_input_vtx(void *draw_ptr, const nb_mshquad_t *const quad,
			   const camera_t *cam)
{
	double r = 3.0;
	nb_drawing_set_source_rgb(draw_ptr, 0.9, 0.2, 0.4);
	for (uint32_t i = 0; i < quad->N_vtx; i++) {
		uint32_t ni = quad->vtx[i];
		if (ni < quad->N_nod) {
			double x = quad->nod[ni * 2];
			double y = quad->nod[ni*2+1];
			nb_drawing_move_to(draw_ptr, cam, x + r, y);
			nb_drawing_arc(draw_ptr, cam, x, y,
				       r, 0, 2 * NB_MATH_PI);
			nb_drawing_fill(draw_ptr);
		}
	}
}
