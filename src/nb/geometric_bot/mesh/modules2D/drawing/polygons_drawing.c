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
			  const void *const poly);
static void draw_input_sgm(nb_graphics_context_t *g,
			   const void *const poly);

void nb_mshpoly_draw(const void *const poly,
		     const char* filename, int width, int height)
{
	nb_graphics_export(filename, width, height, draw_mesh, poly);
}


static void draw_mesh(nb_graphics_context_t *g, int width, int height,
		      const void *const poly)
{
	double box[4];
	nb_mshpoly_get_enveloping_box(poly, box);

	nb_graphics_enable_camera(g);
	nb_graphics_camera_t* cam = nb_graphics_get_camera(g);
	nb_graphics_cam_fit_box(cam, box, width, height);

	draw_polygons(g, poly);
	draw_input_sgm(g, poly);
}

static void draw_polygons(nb_graphics_context_t *g,
			  const void *const poly)
{
	nb_graphics_set_line_width(g, 0.5);
	uint32_t N_elems = nb_mshpoly_get_N_elems(poly);
	for (uint32_t i = 0; i < N_elems; i++) {
		uint32_t id = nb_mshpoly_elem_get_adj(poly, i, 0);
		double x = nb_mshpoly_get_x_node(poly, id);
		double y = nb_mshpoly_get_y_node(poly, id);
		nb_graphics_move_to(g, x, y);
		uint16_t N_adj = nb_mshpoly_elem_get_N_adj(poly, i);
		for (uint16_t j = 1; j < N_adj; j++) {
			id = nb_mshpoly_elem_get_adj(poly, i, j);
			x = nb_mshpoly_get_x_node(poly, id);
			y = nb_mshpoly_get_y_node(poly, id);
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
			   const void *const poly)
{
	nb_graphics_set_line_width(g, 1.0);
	nb_graphics_set_source(g, NB_AQUAMARIN);
	uint32_t N_sgm = nb_mshpoly_get_N_insgm(poly);
	for (uint32_t i = 0; i < N_sgm; i++) {
		uint32_t N_nod_x_sgm = nb_mshpoly_get_N_nodes_x_insgm(poly, i);
		if (0 < N_nod_x_sgm) {
			uint32_t nj = nb_mshpoly_get_node_x_insgm(poly, i, 0);
			nb_graphics_move_to(g,
					    nb_mshpoly_get_x_node(poly, nj),
					    nb_mshpoly_get_y_node(poly, nj));
			for (uint32_t j = 1; j < N_nod_x_sgm; j++) {
				nj = nb_mshpoly_get_node_x_insgm(poly, i, j);
				nb_graphics_line_to(g,
						    nb_mshpoly_get_x_node(poly, nj),
						    nb_mshpoly_get_y_node(poly, nj));
			}
			nb_graphics_stroke(g);
		}
	}

}
