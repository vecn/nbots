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

static void draw_fem(nb_graphics_context_t *g, int width, int height,
		     const void *const data);
static void set_source_trg(nb_graphics_context_t *g,
			   const double v1[2],
			   const double v2[2],
			   const double v3[2],
			   const nb_graphics_palette_t *palette,
			   const double *results,
			   uint32_t min_id, uint32_t max_id,
			   uint32_t id1, uint32_t id2, uint32_t id3);

void nb_fem_save(const vcn_msh3trg_t *const msh3trg,
		 const double *results,
		 const char* filename,
		 int width, int height)
{
	fem_data_t data;
	data.mesh = msh3trg;
	data.results = results;
	nb_graphics_export(filename, width, height, draw_fem,
			   &data);
}

static void draw_fem(nb_graphics_context_t *g, int width, int height,
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
	
	nb_graphics_enable_camera(g);
	nb_graphics_camera_t* cam = nb_graphics_get_camera(g);
	nb_graphics_cam_fit_box(cam, box, width, height);

	nb_graphics_palette_t *palette = 
		nb_graphics_palette_create_preset(NB_RAINBOW);

	uint32_t min_id;
	uint32_t max_id;
	vcn_array_get_min_max_ids(results, msh3trg->N_vertices, sizeof(*results),
				  vcn_compare_double, &min_id, &max_id);

	/* Draw triangles */
	nb_graphics_set_line_width(g, 0.5);
	for (uint32_t i = 0; i < msh3trg->N_triangles; i++) {
		uint32_t n1 = msh3trg->vertices_forming_triangles[i * 3];
		uint32_t n2 = msh3trg->vertices_forming_triangles[i*3+1];
		uint32_t n3 = msh3trg->vertices_forming_triangles[i*3+2];
		double *v1 = &(msh3trg->vertices[n1*2]);
		double *v2 = &(msh3trg->vertices[n2*2]);
		double *v3 = &(msh3trg->vertices[n3*2]);
		nb_graphics_move_to(g, v1[0], v1[1]);
		nb_graphics_line_to(g, v2[0], v2[1]);
		nb_graphics_line_to(g, v3[0], v3[1]);
		nb_graphics_close_path(g);

		set_source_trg(g, v1, v2, v3, palette, results,
			       min_id, max_id, n1, n2, n3);
		nb_graphics_fill_preserve(g);

		nb_graphics_set_source(g, NB_BLACK);
		nb_graphics_stroke(g);
	}

	/* Draw input segments */
	nb_graphics_set_source_rgb(g, 0.9, 0.0, 0.8);

	nb_graphics_set_line_width(g, 1.0);
	for (uint32_t i = 0; i < msh3trg->N_input_segments; i++) {
		if (msh3trg->N_subsgm_x_inputsgm[i] == 0)
			continue;

		uint32_t n1 = msh3trg->meshvtx_x_inputsgm[i][0];
		double *v1 = &(msh3trg->vertices[n1*2]);
		nb_graphics_move_to(g, v1[0], v1[1]);
		for (uint32_t j = 0; j < msh3trg->N_subsgm_x_inputsgm[i]; j++) {
			uint32_t n2 = msh3trg->meshvtx_x_inputsgm[i][j+1];
			double *v2 = &(msh3trg->vertices[n2*2]);
			nb_graphics_line_to(g, v2[0], v2[1]);
		}
		nb_graphics_stroke(g);
	}
}

static void set_source_trg(nb_graphics_context_t *g,
			   const double v1[2],
			   const double v2[2],
			   const double v3[2],
			   const nb_graphics_palette_t *palette,
			   const double *results,
			   uint32_t min_id, uint32_t max_id,
			   uint32_t id1, uint32_t id2, uint32_t id3)
{
	double results_range = results[max_id] - results[min_id];
	uint8_t c1[4], c2[4], c3[4];
	double val = (results[id1] - results[min_id]) / results_range;
	nb_graphics_palette_get_rgba(palette, val, c1);

	val = (results[id2] - results[min_id]) / results_range;
	nb_graphics_palette_get_rgba(palette, val, c2);

	val = (results[id3] - results[min_id]) / results_range;
	nb_graphics_palette_get_rgba(palette, val, c3);

	nb_graphics_set_source_trg(g, v1[0], v1[1], v2[0],
				   v2[1], v3[0], v3[1],
				   c1, c2, c3);
}
