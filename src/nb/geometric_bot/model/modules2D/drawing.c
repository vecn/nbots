#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>

#include "nb/math_bot.h"
#include "nb/container_bot.h"
#include "nb/graphics_bot.h"
#include "nb/geometric_bot.h"

static void draw_edges(void *draw_ptr, const camera_t *cam,
		       const vcn_model_t *const model);
static void draw_vertices(void *draw_ptr, const camera_t *cam,
			  const vcn_model_t *const model);
static void draw_holes(void *draw_ptr, const camera_t *cam,
		       const vcn_model_t *const model);
static void draw_model(void *draw_ptr, int width, int height,
		       const void *const model_ptr);

static void draw_edges(void *draw_ptr, const camera_t *cam,
		       const vcn_model_t *const model)
{
	nb_drawing_set_source_rgb(draw_ptr, 1.0, 0.0, 0.8);
	nb_drawing_set_line_width(draw_ptr, 1.0);
	for (uint32_t i = 0; i < model->M; i++) {
		uint32_t id1 = model->edge[i * 2];
		uint32_t id2 = model->edge[i*2+1];
		nb_drawing_move_to(draw_ptr, cam,
				   model->vertex[id1 * 2],
				   model->vertex[id1*2+1]);
		nb_drawing_line_to(draw_ptr, cam,
				   model->vertex[id2 * 2],
				   model->vertex[id2*2+1]);
		nb_drawing_stroke(draw_ptr);
	}
}

static void draw_vertices(void *draw_ptr, const camera_t *cam,
			  const vcn_model_t *const model)
{
	for (uint32_t i = 0; i < model->N; i++) {
		nb_drawing_set_source_rgb(draw_ptr, 0.0, 1.0, 0.4);

		nb_drawing_move_to(draw_ptr, cam,
				   model->vertex[i * 2],
				   model->vertex[i*2+1]);
		nb_drawing_set_circle(draw_ptr, cam,
				      model->vertex[i * 2],
				      model->vertex[i*2+1],
				      3.0, true);
		nb_drawing_fill(draw_ptr);

		/* Draw labels */
		nb_drawing_set_source_rgb(draw_ptr, 0.0, 0.0, 0.0);

		/* Show id */
		char str_id[5];
		sprintf(str_id, "%i", i);
		nb_drawing_set_font_type(draw_ptr, "Sans");
		nb_drawing_set_font_size(draw_ptr, 9);
		nb_drawing_move_to(draw_ptr, cam,
				   model->vertex[i * 2],
				   model->vertex[i*2+1]);
		nb_drawing_show_text(draw_ptr, str_id);
	}
}

static void draw_holes(void *draw_ptr, const camera_t *cam,
		       const vcn_model_t *const model)
{
	nb_drawing_set_source_rgb(draw_ptr, 0.0, 0.0, 1.0);
	for (uint32_t i = 0; i < model->H; i++) {
		nb_drawing_move_to(draw_ptr, cam,
				   model->holes[i * 2],
				   model->holes[i*2+1]);
		nb_drawing_set_circle(draw_ptr, cam,
				      model->holes[i * 2],
				      model->holes[i*2+1], 
				      3.0, true);
		nb_drawing_fill(draw_ptr);
	}
}

static void draw_model(void *draw_ptr, int width, int height,
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

	camera_t cam;
	nb_drawing_utils_set_center_and_zoom(&cam, box, width, height);

	draw_edges(draw_ptr, &cam, model);
	draw_vertices(draw_ptr, &cam, model);
	draw_holes(draw_ptr, &cam, model);
	
EXIT:
	return;
}

void vcn_model_export_png(const vcn_model_t *const model,
			  const char* filename,
			  int width, int height)
{
	nb_drawing_export_png(filename, width, height, draw_model, model);
}

void vcn_model_export_eps(const vcn_model_t *const model,
			  const char* filename,
			  int width, int height) 
{
	nb_drawing_export_eps(filename, width, height, draw_model, model);
}
