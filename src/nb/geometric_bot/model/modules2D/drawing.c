#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>

#include "nb/math_bot.h"
#include "nb/container_bot.h"
#include "nb/graphics_bot.h"
#include "nb/geometric_bot.h"

static void draw_edges(nb_graphics_context_t *g,
		       const vcn_model_t *const model);
static void draw_vertices(nb_graphics_context_t *g,
			  const vcn_model_t *const model);
static void draw_holes(nb_graphics_context_t *g,
		       const vcn_model_t *const model);
static void draw_model(nb_graphics_context_t *g, int width, int height,
		       const void *const model_ptr);

static void draw_edges(nb_graphics_context_t *g,
		       const vcn_model_t *const model)
{
	nb_graphics_set_source(g, NB_CHARTREUSE);
	nb_graphics_set_line_width(g, 1.0);
	for (uint32_t i = 0; i < model->M; i++) {
		uint32_t id1 = model->edge[i * 2];
		uint32_t id2 = model->edge[i*2+1];
		nb_graphics_move_to(g,
				    model->vertex[id1 * 2],
				    model->vertex[id1*2+1]);
		nb_graphics_line_to(g,
				    model->vertex[id2 * 2],
				    model->vertex[id2*2+1]);
		nb_graphics_stroke(g);
	}
}

static void draw_vertices(nb_graphics_context_t *g,
			  const vcn_model_t *const model)
{
	for (uint32_t i = 0; i < model->N; i++) {
		nb_graphics_set_source(g, NB_AQUAMARIN);

		nb_graphics_move_to(g,
				    model->vertex[i * 2],
				    model->vertex[i*2+1]);
		nb_graphics_set_point(g,
				      model->vertex[i * 2],
				      model->vertex[i*2+1],
				      6.0);
		nb_graphics_fill(g);

		/* Draw id-labels */
		nb_graphics_set_source(g, NB_BLACK);

		char str_id[5];
		sprintf(str_id, "%i", i);
		nb_graphics_set_font_type(g, "Sans");
		nb_graphics_set_font_size(g, 9);
		nb_graphics_show_text(g, 
				      model->vertex[i * 2],
				      model->vertex[i*2+1],
				      str_id);
	}
}

static void draw_holes(nb_graphics_context_t *g,
		       const vcn_model_t *const model)
{
	nb_graphics_set_source(g, NB_BLUE);
	for (uint32_t i = 0; i < model->H; i++) {
		nb_graphics_move_to(g,
				    model->holes[i * 2],
				    model->holes[i*2+1]);
		nb_graphics_set_point(g,
				      model->holes[i * 2],
				      model->holes[i*2+1], 
				      6.0);
		nb_graphics_fill(g);
	}
}

static void draw_model(nb_graphics_context_t *g, int width, int height,
		       const void *const model_ptr)
{
	const vcn_model_t *const model = model_ptr;
	if (0 == model->N) 
		goto EXIT;

	/* Compute center and zoom */
	double box[4];
	vcn_utils2D_get_enveloping_box(model->N, model->vertex,
				       2 * sizeof(*(model->vertex)),
				       vcn_utils2D_get_x_from_darray,
				       vcn_utils2D_get_y_from_darray,
				       box);

	nb_graphics_enable_camera(g);
	nb_graphics_camera_t* cam = nb_graphics_get_camera(g);
	nb_graphics_cam_fit_box(cam, box, width, height);

	draw_edges(g, model);
	draw_vertices(g, model);
	draw_holes(g, model);
	
EXIT:
	return;
}

void vcn_model_draw(const vcn_model_t *const model,
		    const char* filename,
		    int width, int height)
{
	nb_graphics_export(filename, width, height, draw_model, model);
}
