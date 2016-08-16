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
		      const void *const quad);
static void draw_quads(nb_graphics_context_t *g,
		       const void *const quad);
static void draw_input_sgm(nb_graphics_context_t *g,
			   const void *const quad);
static void draw_input_vtx(nb_graphics_context_t *g,
			   const void *const quad);

void nb_mshquad_draw(const void *const quad,
		     const char* filename, int width, int height)
{
	nb_graphics_export(filename, width, height, draw_mesh, quad);
}

static void draw_mesh(nb_graphics_context_t *g, int width, int height,
		      const void *const quad)
{

	double box[4];
	nb_mshquad_get_enveloping_box(quad, box);

	nb_graphics_enable_camera(g);
	nb_graphics_camera_t* cam = nb_graphics_get_camera(g);
	nb_graphics_cam_fit_box(cam, box, width, height);

	draw_quads(g, quad);
	draw_input_sgm(g, quad);
	draw_input_vtx(g, quad);
}

static void draw_quads(nb_graphics_context_t *g,
		       const void *const quad)
{
	nb_graphics_set_line_width(g, 0.5);
	uint32_t N_elems = nb_mshquad_get_N_elems(quad);
	for (uint32_t i = 0; i < N_elems; i++) {
		uint32_t n1 = nb_mshquad_elem_get_adj(quad, i, 0);
		uint32_t n2 = nb_mshquad_elem_get_adj(quad, i, 1);
		uint32_t n3 = nb_mshquad_elem_get_adj(quad, i, 2);
		uint32_t n4 = nb_mshquad_elem_get_adj(quad, i, 3);
		nb_graphics_move_to(g,
				    nb_mshquad_get_x_node(quad, n1),
				    nb_mshquad_get_y_node(quad, n1));
		nb_graphics_line_to(g,
				    nb_mshquad_get_x_node(quad, n2),
				    nb_mshquad_get_y_node(quad, n2));
		nb_graphics_line_to(g,
				    nb_mshquad_get_x_node(quad, n3),
				    nb_mshquad_get_y_node(quad, n3));
		if (n4 < nb_mshquad_get_N_nodes(quad))
			nb_graphics_line_to(g,
					    nb_mshquad_get_x_node(quad, n4),
					    nb_mshquad_get_y_node(quad, n4));
		nb_graphics_close_path(g);

		nb_graphics_set_source_rgba(g, 25, 75, 255, 128);
		nb_graphics_fill_preserve(g);

		nb_graphics_set_source(g, NB_BLUE);
		nb_graphics_stroke(g);
	}
}

static void draw_input_sgm(nb_graphics_context_t *g, const void *const quad)
{
	nb_graphics_set_line_width(g, 1.0);
	nb_graphics_set_source(g, NB_AQUAMARIN);
	uint32_t N_sgm = nb_mshquad_get_N_insgm(quad);
	for (uint32_t i = 0; i < N_sgm; i++) {
		uint32_t N_nod_x_sgm = nb_mshquad_get_N_nodes_x_insgm(quad, i);
		if (0 < N_nod_x_sgm) {
			uint32_t nj = nb_mshquad_get_node_x_insgm(quad, i, 0);
			nb_graphics_move_to(g,
					    nb_mshquad_get_x_node(quad, nj),
					    nb_mshquad_get_y_node(quad, nj));
			for (uint32_t j = 1; j < N_nod_x_sgm; j++) {
				nj = nb_mshquad_get_node_x_insgm(quad, i, j);
				double x = nb_mshquad_get_x_node(quad, nj);
				double y = nb_mshquad_get_y_node(quad, nj);
				nb_graphics_line_to(g, x, y);
			}
			nb_graphics_stroke(g);
		}
	}
}

static void draw_input_vtx(nb_graphics_context_t *g,
			   const void *const quad)
{
	nb_graphics_set_source(g, NB_CHARTREUSE);
	uint32_t N_vtx = nb_mshquad_get_N_invtx(quad);
	for (uint32_t i = 0; i < N_vtx; i++) {
		uint32_t ni = nb_mshquad_get_invtx(quad, i);
		if (ni < nb_mshquad_get_N_nodes(quad)) {
			double x = nb_mshquad_get_x_node(quad, ni);
			double y = nb_mshquad_get_y_node(quad, ni);
			nb_graphics_set_point(g, x, y, 6.0);
			nb_graphics_fill(g);
		}
	}
}
