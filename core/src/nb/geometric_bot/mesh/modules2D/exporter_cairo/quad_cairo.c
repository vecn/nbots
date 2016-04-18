#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <cairo.h>
#include <cairo-ps.h>

#include "nb/math_bot.h"
#include "nb/container_bot.h"
#include "nb/geometric_bot.h"
#include "nb/geometric_bot-cairo.h"

#include "../../mesh2D_structs.h"
#include "drawing_utils.h"
#include "quad_cairo.h"

static void draw_mesh(cairo_t *cr, int width, int height,
		      const nb_mshquad_t *const quad);
static void set_camera(camera_t *cam, const nb_mshquad_t *const quad,
		       int width, int height);

void nb_mshquad_export_png(const nb_mshquad_t *const quad,
			   const char* filename, int width, int height)
{
	nb_cairo_drawing_export_png(filename, width, height, draw_mesh, quad);
}

void nb_mshquad_export_eps(const nb_mshquad_t *const quad,
			   const char* filename, int width, int height)
{
	nb_cairo_drawing_export_eps(filename, width, height, draw_mesh, quad);
}

static void draw_mesh(cairo_t *cr, int width, int height,
		      const nb_mshquad_t *const quad)
{
	camera_t cam;
	set_camera(&cam, quad, width, height);

	cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
	cairo_paint(cr);

	double rgba_bg[4] = {0.1, 0.1, 1.0, 0.5};

	cairo_set_line_width(cr, 0.5);
	draw_quads(cr, quad, width, height, cam, 
		       _NB_COLOR_BLUE);

	double rgb_segments[3] = {1.0, 0.0, 0.8};
	draw_input_segments(cr, mesh->N_input_sgm, mesh->input_sgm,
			    width, height, cam, rgb_segments);

	draw_input_vertices(cr, mesh->N_input_vtx, mesh->input_vtx,
			    width, height, cam, false,
			    _NB_COLOR_BLUE, _NB_COLOR_BLACK);
}

static void set_camera(camera_t *cam, const nb_mshquad_t *const quad,
		       int width, int height)
{

}
