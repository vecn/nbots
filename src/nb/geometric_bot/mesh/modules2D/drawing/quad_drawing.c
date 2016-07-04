#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>

#include "nb/math_bot.h"
#include "nb/container_bot.h"
#include "nb/geometric_bot.h"
#include "nb/graphics_bot.h"

static void draw_mesh(nb_graphics_context_t *g, int width, int height,
		      const void *const quad_ptr);
static void draw_quads(nb_graphics_context_t *g,
		       const nb_mshquad_t *const quad);
static void draw_input_sgm(nb_graphics_context_t *g,
			   const nb_mshquad_t *const quad);
static void draw_input_vtx(nb_graphics_context_t *g,
			   const nb_mshquad_t *const quad);

void nb_mshquad_draw(const nb_mshquad_t *const quad,
		     const char* filename, int width, int height)
{
	nb_graphics_export(filename, width, height, draw_mesh, quad);
}

static void draw_mesh(nb_graphics_context_t *g, int width, int height,
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
	nb_graphics_enable_camera(g);
	nb_graphics_camera_t* cam = nb_graphics_get_camera(g);
	nb_graphics_cam_fit_box(cam, box, width, height);

	draw_quads(g, quad);
	draw_input_sgm(g, quad);
	draw_input_vtx(g, quad);
}

static void draw_quads(nb_graphics_context_t *g,
		       const nb_mshquad_t *const quad)
{
	nb_graphics_set_line_width(g, 0.5);
	for (uint32_t i = 0; i < quad->N_elems; i++) {
		uint32_t n1 = quad->adj[i * 4];
		uint32_t n2 = quad->adj[i*4+1];
		uint32_t n3 = quad->adj[i*4+2];
		uint32_t n4 = quad->adj[i*4+3];
		nb_graphics_move_to(g,
				    quad->nod[n1 * 2],
				    quad->nod[n1*2+1]);
		nb_graphics_line_to(g,
				    quad->nod[n2 * 2],
				    quad->nod[n2*2+1]);
		nb_graphics_line_to(g,
				    quad->nod[n3 * 2],
				    quad->nod[n3*2+1]);
		if (n4 < quad->N_nod)
			nb_graphics_line_to(g,
					    quad->nod[n4 * 2],
					    quad->nod[n4*2+1]);
		nb_graphics_close_path(g);

		nb_graphics_set_source_rgba(g, 25, 75, 255, 128);
		nb_graphics_fill_preserve(g);

		nb_graphics_set_source(g, NB_BLUE);
		nb_graphics_stroke(g);
	}
}

static void draw_input_sgm(nb_graphics_context_t *g,
			   const nb_mshquad_t *const quad)
{
	nb_graphics_set_line_width(g, 1.0);
	nb_graphics_set_source(g, NB_AQUAMARIN);
	for (uint32_t i = 0; i < quad->N_sgm; i++) {
		if (0 < quad->N_nod_x_sgm[i]) {
			uint32_t nj = quad->nod_x_sgm[i][0];
			nb_graphics_move_to(g,
					    quad->nod[nj * 2],
					    quad->nod[nj*2+1]);
			for (uint32_t j = 1; j < quad->N_nod_x_sgm[i]; j++) {
				nj = quad->nod_x_sgm[i][j];
				nb_graphics_line_to(g,
						    quad->nod[nj * 2],
						    quad->nod[nj*2+1]);
			}
			nb_graphics_stroke(g);
		}
	}
}

static void draw_input_vtx(nb_graphics_context_t *g,
			   const nb_mshquad_t *const quad)
{
	nb_graphics_set_source(g, NB_CHARTREUSE);
	for (uint32_t i = 0; i < quad->N_vtx; i++) {
		uint32_t ni = quad->vtx[i];
		if (ni < quad->N_nod) {
			double x = quad->nod[ni * 2];
			double y = quad->nod[ni*2+1];
			nb_graphics_set_point(g, x, y, 6.0);
			nb_graphics_fill(g);
		}
	}
}
