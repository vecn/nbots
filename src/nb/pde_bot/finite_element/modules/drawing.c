#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "nb/math_bot.h"
#include "nb/container_bot.h"
#include "nb/geometric_bot.h"
#include "nb/graphics_bot.h"

#include "nb/pde_bot/finite_element/modules/drawing.h"

#define POW2(a) ((a)*(a))

typedef struct {
	const nb_partition_t *part;
	const double *distortion;
	double max_distortion;
	const double *results;
} fem_data_t;

static void draw_fem(nb_graphics_context_t *g, int width, int height,
		     const void *const data);

static double get_distortion_scale(const nb_partition_t *part,
				   const double *distortion,
				   double max_distortion);
static void get_enveloping_box(const nb_partition_t *part,
			       const double *distortion,
			       double dscale, double box[4]);
static void get_distortion(double x[2], uint32_t vtx_id,
			   const nb_partition_t *part,
			   const double *distortion,
			   double dscale);
static void set_source_trg(nb_graphics_context_t *g,
			   const double v1[2],
			   const double v2[2],
			   const double v3[2],
			   const nb_graphics_palette_t *palette,
			   const double *results,
			   uint32_t min_id, uint32_t max_id,
			   uint32_t id1, uint32_t id2, uint32_t id3);

void nb_fem_save(const nb_partition_t *const part,
		 const double *distortion, /* Can be NULL */
		 double max_distortion,
		 const double *results,
		 const char* filename,
		 int width, int height)
{
	fem_data_t data;
	data.part = part;
	data.results = results;
	data.distortion = distortion;
	data.max_distortion = max_distortion;
	nb_graphics_export(filename, width, height, draw_fem,
			   &data);
}

static void draw_fem(nb_graphics_context_t *g, int width, int height,
		     const void *const data)
{
	const fem_data_t *const fem_data = data;
	const nb_partition_t *part = fem_data->part;
	const double *results = fem_data->results;
	const double *distortion = fem_data->distortion;
	double max_distortion = fem_data->max_distortion;

	double dscale = get_distortion_scale(part,
					     distortion,
					     max_distortion);
	
	double box[4];
	get_enveloping_box(part, distortion, dscale, box);
	
	nb_graphics_enable_camera(g);
	nb_graphics_camera_t* cam = nb_graphics_get_camera(g);
	nb_graphics_cam_fit_box(cam, box, width, height);

	nb_graphics_palette_t *palette = 
		nb_graphics_palette_create_preset(NB_RAINBOW);

	uint32_t min_id;
	uint32_t max_id;
	uint32_t N_nod = nb_partition_get_N_nodes(part);
	vcn_array_get_min_max_ids(results, N_nod, sizeof(*results),
				  vcn_compare_double, &min_id, &max_id);

	/* Draw triangles */
	nb_graphics_set_line_width(g, 0.5);
	uint32_t N_elems = nb_partition_get_N_elems(part);
	for (uint32_t i = 0; i < N_elems; i++) {
		uint32_t n1 = nb_partition_elem_get_adj(part, i, 0);
		uint32_t n2 = nb_partition_elem_get_adj(part, i, 1);
		uint32_t n3 = nb_partition_elem_get_adj(part, i, 2);

		double v1[2], v2[2], v3[2];
		get_distortion(v1, n1, part, distortion, dscale);
		get_distortion(v2, n2, part, distortion, dscale);
		get_distortion(v3, n3, part, distortion, dscale);

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
	uint32_t N_sgm = nb_partition_get_N_insgm(part);
	for (uint32_t i = 0; i < N_sgm; i++) {
		uint32_t N_vtx = nb_partition_get_N_nodes_x_insgm(part, i);
		if (0 < N_vtx) {
			uint32_t n1 = nb_partition_get_node_x_insgm(part, i, 0);
			double v1[2];
			get_distortion(v1, n1, part,
				       distortion, dscale);
			nb_graphics_move_to(g, v1[0], v1[1]);
			for (uint32_t j = 0; j < N_vtx; j++) {
				uint32_t n2 = nb_partition_get_node_x_insgm(part, i, j);
				double v2[2];
				get_distortion(v2, n2, part,
					       distortion, dscale);
				nb_graphics_line_to(g, v2[0], v2[1]);
			}
			nb_graphics_stroke(g);
		}
	}

	nb_graphics_palette_draw(g, palette, width - 100, height - 400,
				 20, 300, 1.0f,
				 results[min_id], results[max_id]);
	nb_graphics_palette_destroy(palette);
}

static double get_distortion_scale(const nb_partition_t *part,
				   const double *distortion,
				   double max_distortion)
{
	double dscale;
	if (NULL != distortion) {
		double max_dist = 0.0;
		uint32_t N_sgm = nb_partition_get_N_insgm(part);
		for (uint32_t i = 0; i < N_sgm; i++) {
			uint32_t N_vtx = nb_partition_get_N_nodes_x_insgm(part, i);
			for (uint32_t j = 0; j < N_vtx; j++) {
				uint32_t id = nb_partition_get_node_x_insgm(part, i, j);
				double dist = POW2(distortion[id * 2]) +
					POW2(distortion[id*2+1]);
				if (dist > max_dist)
					max_dist = dist;
			}
		}
		dscale = max_distortion / sqrt(max_dist);
	} else {
		dscale = 0.0;
	}
	return dscale;
}

static void get_enveloping_box(const nb_partition_t *part,
			       const double *distortion,
			       double dscale, double box[4])
{	
	uint32_t N_sgm = nb_partition_get_N_insgm(part);
	for (uint32_t i = 0; i < N_sgm; i++) {
		uint32_t N_vtx = nb_partition_get_N_nodes_x_insgm(part, i);
		for (uint32_t j = 0; j < N_vtx; j++) {
			uint32_t id = nb_partition_get_node_x_insgm(part, i, j);
			double v[2];
			get_distortion(v, id, part, distortion, dscale);
			if (0 == i && 0 == j) {
				box[0] = v[0];
				box[1] = v[1];
				box[2] = v[0];
				box[3] = v[1];
			} else {
				if (v[0] < box[0])
					box[0] = v[0];
				else if (v[0] > box[2])
					box[2] = v[0];

				if (v[1] < box[1])
					box[1] = v[1];
				else if (v[1] > box[3])
					box[3] = v[1];
			}
		}
	}
}

static void get_distortion(double x[2], uint32_t vtx_id,
			   const nb_partition_t *part,
			   const double *distortion,
			   double dscale)
{
	x[0] = nb_partition_get_x_node(part, vtx_id);
	x[1] = nb_partition_get_y_node(part, vtx_id);
	if (NULL != distortion) {
		x[0] += dscale * distortion[vtx_id * 2];
		x[1] += dscale * distortion[vtx_id*2+1];
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
