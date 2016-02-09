#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <cairo.h>
#include <cairo-ps.h>

#include "vcn/math_bot.h"
#include "vcn/container_bot.h"
#include "vcn/geometric_bot.h"
#include "vcn/geometric_bot-cairo.h"

#include "../model2D_struct.h"

static inline double get_center_and_zoom(double box[4], double width,
					 double height, double center[2]);
static void draw_edges(cairo_t *const cr,
		       uint32_t N_vertices,
		       const double *const vertices,
		       uint32_t N_edges,
		       const uint32_t *const edges,
		       int width, int height,
		       double center[2],
		       double zoom);
static void draw_vertices(cairo_t *const cr,
			  uint32_t N_vertices,
			  const double *const vertices,
			  int width, int height,
			  double center[2], double zoom,
			  bool include_numbering,
			  double rgb_bg[3], double rgb_fg[3]);
static void draw_holes(cairo_t *const cr,
		       uint32_t N_holes,
		       const double *const holes,
		       int width, int height,
		       double center[2], double zoom);
static void model_draw_with_cairo(const vcn_model_t *const model,
				  cairo_t *cr,
				  int width, int height,
				  bool include_numbering);


static inline double get_center_and_zoom(double box[4], double width,
					 double height, double center[2])
{
	center[0] = (box[0] + box[2]) / 2.0;
	center[1] = (box[1] + box[3]) / 2.0;
	double zoom = width / (box[2] - box[0]);
	if (zoom > height / (box[3] - box[1]))
		zoom = height / (box[3] - box[1]);
	return zoom * 0.9;
}

static void draw_edges(cairo_t *const restrict cr,
		       uint32_t N_vertices,
		       const double *const restrict vertices,
		       uint32_t N_edges,
		       const uint32_t *const restrict edges,
		       int width, int height,
		       double center[2],
		       double zoom)
{
	for (uint32_t i = 0; i < N_edges; i++) {
		uint32_t id1 = edges[i * 2];
		uint32_t id2 = edges[i*2+1];
		cairo_move_to(cr,
			      zoom * vertices[id1 * 2] -
			      zoom * center[0] + width/2.0,
			      -zoom * vertices[id1*2+1] +
			      zoom * center[1] + height/2.0);
		cairo_line_to(cr,
			      zoom * vertices[id2 * 2] -
			      zoom * center[0] + width/2.0,
			      -zoom * vertices[id2*2+1] +
			      zoom * center[1] + height/2.0);
		cairo_stroke(cr);
	}
}

static void draw_vertices(cairo_t *const restrict cr,
			  uint32_t N_vertices,
			  const double *const restrict vertices,
			  int width, int height,
			  double center[2], double zoom,
			  bool include_numbering,
			  double rgb_bg[3], double rgb_fg[3])
{
	for (uint32_t i = 0; i < N_vertices; i++) {
		cairo_set_source_rgb(cr, rgb_bg[0], rgb_bg[1], rgb_bg[2]);

		cairo_move_to(cr,
			      zoom * vertices[i * 2] - zoom * center[0] +
			      width / 2.0 + 3.0,
			      -zoom * vertices[i*2+1] + zoom * center[1] +
			      height / 2.0);
		cairo_arc(cr,
			  zoom * vertices[i * 2] - zoom * center[0] +
			  width / 2.0,
			  -zoom * vertices[i*2+1] + zoom * center[1] +
			  height / 2.0, 3.0, 0.0, 2.0 * VCN_MATH_PI);
		cairo_fill(cr);

		/* Draw labels */
		if (include_numbering) {
			cairo_set_source_rgb(cr, rgb_fg[0], rgb_fg[1], rgb_fg[2]);

			/* Show id */
			char str_id[5];
			sprintf(str_id, "%i", i);
			cairo_select_font_face(cr, "Sans",
					       CAIRO_FONT_SLANT_NORMAL,
					       CAIRO_FONT_WEIGHT_NORMAL);
			cairo_set_font_size(cr, 9);
			cairo_move_to(cr, 
				      zoom * vertices[i * 2] - 
				      zoom * center[0] + width/2.0 + 4.0,
				      -zoom * vertices[i*2+1] + 
				      zoom * center[1] + height/2.0 - 1.0);
			cairo_show_text(cr, str_id);
		}
	}
}

static void draw_holes(cairo_t *const restrict cr,
		       uint32_t N_holes,
		       const double *const restrict holes,
		       int width, int height,
		       double center[2], double zoom)
{
	for (uint32_t i = 0; i < N_holes; i++) {
		cairo_move_to(cr,
			      zoom * holes[i * 2] - zoom * center[0] +
			      width / 2.0 + 3.0,
			      -zoom * holes[i*2+1] + zoom * center[1] +
			      height / 2.0);
		cairo_arc(cr,
			  zoom * holes[i * 2] - zoom * center[0] +
			  width / 2.0,
			  -zoom * holes[i*2+1] + zoom * center[1] +
			  height / 2.0, 
			  3.0, 0.0, 2.0 * VCN_MATH_PI);
		cairo_fill(cr);
	}
}

static void model_draw_with_cairo(const vcn_model_t *const model,
				  cairo_t* cr,
				  int width, int height,
				  bool include_numbering) 
{
	if (model->N == 0)
		return;
  
	/* Compute center and zoom */
	double box[4];
	vcn_utils2D_get_enveloping_box(model->N, model->vertex,
				       2 * sizeof(*(model->vertex)),
				       vcn_utils2D_get_x_from_darray,
				       vcn_utils2D_get_y_from_darray,
				       box);

	double center[2];
	double zoom = get_center_and_zoom(box, width, height, center);

	/* Draw background */
	cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
	cairo_paint(cr);

	/* Draw input segments */
	cairo_set_source_rgb(cr, 1.0, 0.0, 0.8);
	cairo_set_line_width(cr, 0.75);
	draw_edges(cr, model->N, model->vertex, model->M, model->edge,
		   width, height, center, zoom);

	/* Draw input vertices */
	double rgb_bg[3] = {1.0, 0.0, 0.0};
	double rgb_fg[3] = {0.0, 0.0, 0.0};
	cairo_set_line_width(cr, 0.75);
	draw_vertices(cr, model->N, model->vertex, width, height,
		      center, zoom, include_numbering, rgb_bg, rgb_fg);

	/* Draw holes */
	cairo_set_source_rgb(cr, 0.0, 0.0, 1.0);
	draw_holes(cr, model->H, model->holes, width, height, center, zoom);
}

void vcn_model_export_png(const vcn_model_t *const model,
			  const char* filename,
			  int width, int height,
			  bool include_numbering)
{
  
	/* Create drawable surface and cairo context */
	cairo_surface_t* surface =
		cairo_image_surface_create(CAIRO_FORMAT_ARGB32, 
					   width, height);
	cairo_t *const restrict cr = cairo_create(surface);

	/* Draw with cairo */
	model_draw_with_cairo(model, cr, width, height, include_numbering);

	/* Write file */
	cairo_surface_write_to_png(surface, filename);

	/* Free cairo structures */
	cairo_destroy(cr);
	cairo_surface_destroy(surface);
}

void vcn_model_export_eps(const vcn_model_t *const model,
		      const char* filename,
		      int width, int height,
		      bool include_numbering) 
{
	/* Create drawable surface and cairo context */
	cairo_surface_t* surface = 
		cairo_ps_surface_create (filename, width, height);
	cairo_ps_surface_set_eps(surface, 1 /* TRUE from cairo_bool_t*/);

	cairo_t* cr = cairo_create(surface);

	/* Initialize Post script */
	cairo_ps_surface_dsc_begin_page_setup(surface);

	if (height < width)
		cairo_ps_surface_dsc_comment(surface,
					     "%%PageOrientation: Portrait");
	else
		cairo_ps_surface_dsc_comment(surface,
					     "%%PageOrientation: Landscape");

	/* Draw with cairo */
	model_draw_with_cairo(model, cr, width, height, include_numbering);

	/* Show post-script */
	cairo_surface_show_page(surface);

	/* Free cairo structures */
	cairo_destroy(cr);
	cairo_surface_finish(surface);
	cairo_surface_destroy(surface);
}
