#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>

#include "nb/math_bot.h"
#include "nb/container_bot.h"
#include "nb/geometric_bot.h"
#include "nb/graphics_bot.h"

#include "nb/pde_bot/finite_element/modules/drawing.h"

typedef struct {
	const vcn_msh3trg_t *mesh;
	const double *results;
} fem_data_t;

static void draw_fem(void *draw_ptr, int width, int height,
		     const void *const data);
static void set_pattern_patch(nb_pattern_t *pat, const camera_t *const cam,
			      const double v1[2], const double v2[2],
			      const double v3[2]);
static void set_pattern_color(nb_pattern_t *pat,
			      const vcn_palette_t *palette,
			      const double *results,
			      uint32_t min_id, uint32_t max_id,
			      uint32_t id1, uint32_t id2, uint32_t id3);

void nb_fem_save_png(const vcn_msh3trg_t *const msh3trg,
		     const double *results,
		     const char* filename,
		     int width, int height)
{
	fem_data_t data;
	data.mesh = msh3trg;
	data.results = results;
	nb_drawing_export_png(filename, width, height, draw_fem,
			      &data);
}

static void draw_fem(void *draw_ptr, int width, int height,
		     const void *const data)
{
	const fem_data_t *const fem_data = data;
	const vcn_msh3trg_t *msh3trg = fem_data->mesh;
	const double *results = fem_data->results;

	/* Compute cam->center and cam->zoom */
	double box[4];
	vcn_utils2D_get_enveloping_box_from_subset
		(msh3trg->N_input_vertices,
		 msh3trg->input_vertices,
		 msh3trg->vertices,
		 2 * sizeof(*(msh3trg->vertices)),
		 vcn_utils2D_get_x_from_darray,
		 vcn_utils2D_get_y_from_darray,
		 box);

	camera_t cam;
	nb_drawing_utils_set_center_and_zoom(&cam, box, width, height);

	vcn_palette_t *palette = vcn_palette_create_preset(NB_RAINBOW);

	uint32_t min_id;
	uint32_t max_id;
	vcn_array_get_min_max_ids(results, msh3trg->N_vertices, sizeof(*results),
				  vcn_compare_double, &min_id, &max_id);

	/* Draw triangles */
	nb_drawing_set_line_width(draw_ptr, 0.5);
	for (uint32_t i = 0; i < msh3trg->N_triangles; i++) {
		uint32_t n1 = msh3trg->vertices_forming_triangles[i * 3];
		uint32_t n2 = msh3trg->vertices_forming_triangles[i*3+1];
		uint32_t n3 = msh3trg->vertices_forming_triangles[i*3+2];
		double *v1 = &(msh3trg->vertices[n1*2]);
		double *v2 = &(msh3trg->vertices[n2*2]);
		double *v3 = &(msh3trg->vertices[n3*2]);
		nb_drawing_move_to(draw_ptr, &cam, v1[0], v1[1]);
		nb_drawing_line_to(draw_ptr, &cam, v2[0], v2[1]);
		nb_drawing_line_to(draw_ptr, &cam, v3[0], v3[1]);
		nb_drawing_close_path(draw_ptr);

		nb_pattern_t *pat = nb_pattern_create();
		nb_pattern_begin_patch(pat);
		set_pattern_patch(pat, &cam, v1, v2, v3);
		set_pattern_color(pat, palette, results,
				  min_id, max_id, n1, n2, n3);
		nb_pattern_end_patch(pat);
		nb_drawing_set_source(draw_ptr, pat);
		nb_drawing_fill_preserve(draw_ptr);
		nb_pattern_destroy(pat);

		nb_drawing_set_source_rgb(draw_ptr, 0.0, 0.0, 0.0);
		nb_drawing_stroke(draw_ptr);
	}

	/* Draw input segments */
	nb_drawing_set_source_rgb(draw_ptr, 0.9, 0.0, 0.8);

	nb_drawing_set_line_width(draw_ptr, 1.0);
	for (uint32_t i = 0; i < msh3trg->N_input_segments; i++) {
		if (msh3trg->N_subsgm_x_inputsgm[i] == 0)
			continue;

		uint32_t n1 = msh3trg->meshvtx_x_inputsgm[i][0];
		double *v1 = &(msh3trg->vertices[n1*2]);
		nb_drawing_move_to(draw_ptr, &cam, v1[0], v1[1]);
		for (uint32_t j = 0; j < msh3trg->N_subsgm_x_inputsgm[i]; j++) {
			uint32_t n2 = msh3trg->meshvtx_x_inputsgm[i][j+1];
			double *v2 = &(msh3trg->vertices[n2*2]);
			nb_drawing_line_to(draw_ptr, &cam, v2[0], v2[1]);
		}
		nb_drawing_stroke(draw_ptr);
	}
}

static void set_pattern_patch(nb_pattern_t *pat, const camera_t *const cam,
			      const double v1[2], const double v2[2],
			      const double v3[2])
{
	nb_pattern_move_to(pat, cam, v1[0], v1[1]);
	nb_pattern_line_to(pat, cam, v2[0], v2[1]);
	nb_pattern_line_to(pat, cam, v3[0], v3[1]);
}

static void set_pattern_color(nb_pattern_t *pat,
			      const vcn_palette_t *palette,
			      const double *results,
			      uint32_t min_id, uint32_t max_id,
			      uint32_t id1, uint32_t id2, uint32_t id3)
{
	double results_range = results[max_id] - results[min_id];
	uint8_t rgb[3];
	double val = (results[id1] - results[min_id]) / results_range;
	vcn_palette_get_colour(palette, val, rgb);
	nb_pattern_set_corner_color_rgb(pat, 0, 
					rgb[0]/255.0f, 
					rgb[1]/255.0f, 
					rgb[2]/255.0f);
	val = (results[id2] - results[min_id]) / results_range;
	vcn_palette_get_colour(palette, val, rgb);
	nb_pattern_set_corner_color_rgb(pat, 1, 
					rgb[0]/255.0f,
					rgb[1]/255.0f, 
					rgb[2]/255.0f);
	val = (results[id3] - results[min_id]) / results_range;
	vcn_palette_get_colour(palette, val, rgb);
	nb_pattern_set_corner_color_rgb(pat, 2, 
					rgb[0]/255.0f, 
					rgb[1]/255.0f,
					rgb[2]/255.0f);
}
