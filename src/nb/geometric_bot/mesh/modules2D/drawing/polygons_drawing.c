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
		      const void *const poly_ptr);
static void draw_polygons(nb_graphics_context_t *g,
			  const nb_mshpoly_t *const poly);
static void draw_input_sgm(nb_graphics_context_t *g,
			   const nb_mshpoly_t *const poly);

void nb_mshpoly_draw(const nb_mshpoly_t *const poly,
		     const char* filename, int width, int height)
{
	nb_graphics_export(filename, width, height, draw_mesh, poly);
}


static void draw_mesh(nb_graphics_context_t *g, int width, int height,
		      const void *const poly_ptr)
{
	const nb_mshpoly_t *const poly = poly_ptr;

	double box[4];
	vcn_utils2D_get_enveloping_box(poly->N_nod, poly->nod,
				       2 * sizeof(*(poly->nod)),
				       vcn_utils2D_get_x_from_darray,
				       vcn_utils2D_get_y_from_darray,
				       box);
	nb_graphics_enable_camera(g);
	nb_graphics_camera_t* cam = nb_graphics_get_camera(g);
	nb_graphics_cam_fit_box(cam, box, width, height);

	draw_polygons(g, poly);
	draw_input_sgm(g, poly);
}

static void draw_polygons(nb_graphics_context_t *g,
			  const nb_mshpoly_t *const poly)
{
	nb_graphics_set_line_width(g, 0.5);
	for (uint32_t i = 0; i < poly->N_elems; i++) {
		uint32_t id = poly->adj[i][0];
		double x = poly->nod[id * 2];
		double y = poly->nod[id*2+1];
		nb_graphics_move_to(g, x, y);
		for (uint16_t j = 1; j < poly->N_adj[i]; j++) {
			id = poly->adj[i][j];
			x = poly->nod[id * 2];
			y = poly->nod[id*2+1];
			nb_graphics_line_to(g, x, y);			
		}
		nb_graphics_close_path(g);

		nb_graphics_set_source_rgba(g, 25, 75, 255, 128);
		nb_graphics_fill_preserve(g);

		nb_graphics_set_source(g, NB_BLUE);
		nb_graphics_stroke(g);
	}
}

static void draw_input_sgm(nb_graphics_context_t *g,
			   const nb_mshpoly_t *const poly)
{
	nb_graphics_set_line_width(g, 1.0);
	nb_graphics_set_source(g, NB_AQUAMARIN);
	for (uint32_t i = 0; i < poly->N_sgm; i++) {
		if (0 < poly->N_nod_x_sgm[i]) {
			uint32_t nj = poly->nod_x_sgm[i][0];
			nb_graphics_move_to(g,
					    poly->nod[nj * 2],
					    poly->nod[nj*2+1]);
			for (uint32_t j = 1; j < poly->N_nod_x_sgm[i]; j++) {
				nj = poly->nod_x_sgm[i][j];
				nb_graphics_line_to(g,
						    poly->nod[nj * 2],
						    poly->nod[nj*2+1]);
			}
			nb_graphics_stroke(g);
		}
	}

}
