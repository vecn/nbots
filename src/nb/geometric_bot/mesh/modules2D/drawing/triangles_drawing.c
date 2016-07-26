#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>

#include "nb/math_bot.h"
#include "nb/container_bot.h"
#include "nb/geometric_bot.h"
#include "nb/graphics_bot.h"

typedef struct {
	const vcn_msh3trg_t *msh3trg;
	uint32_t k_part;
	const uint32_t *part;
	uint32_t k_to_draw;
	double scale;
} part_data_t;

static void draw_msh3trg(nb_graphics_context_t *g, int width, int height,
			 const void *const msh3trg_ptr);
static void calculate_partition_centers(const vcn_msh3trg_t *const msh3trg,
					const uint32_t *const part,
					uint32_t kpart, double **pcenter);
static void draw_triangle_partition_edge(nb_graphics_context_t *g,
					 const vcn_msh3trg_t *const msh3trg,
					 const uint32_t *const part,
					 uint32_t kpart, double **pcenter,
					 double scale_partitions,
					 int trg_id, int edge_id);
static void draw_triangle_partition_border(nb_graphics_context_t *g,
					   const vcn_msh3trg_t *const msh3trg,
					   const uint32_t *const part,
					   uint32_t kpart, double **pcenter,
					   double scale_partitions,
					   uint32_t trg_id);

static void draw_msh3trg_partition(nb_graphics_context_t *g,
				   int width, int height,
				   const void *const part_data);

static void scale_vtx(double vtx[2], double center[2], double zoom);

void vcn_msh3trg_draw(const vcn_msh3trg_t *const msh3trg,
		      const char* filename, int width, int height)
{
	nb_graphics_export(filename, width, height,
			   draw_msh3trg, msh3trg);
}

static void draw_msh3trg(nb_graphics_context_t *g, int width, int height,
			 const void *const msh3trg_ptr)
{
	const vcn_msh3trg_t *const msh3trg = msh3trg_ptr;
	
	/* Compute cam->center and cam->zoom */
	double box[4];
	vcn_utils2D_get_enveloping_box_from_subset(msh3trg->N_input_vertices,
						   msh3trg->input_vertices,
						   msh3trg->vertices,
						   2 * sizeof(*(msh3trg->vertices)),
						   vcn_utils2D_get_x_from_darray,
						   vcn_utils2D_get_y_from_darray,
						   box);

	nb_graphics_enable_camera(g);
	nb_graphics_camera_t* cam = nb_graphics_get_camera(g);
	nb_graphics_cam_fit_box(cam, box, width, height);

	/* Draw triangles */
	nb_graphics_set_line_width(g, 0.5);
	for (uint32_t i = 0; i < msh3trg->N_triangles; i++) {
		uint32_t n1 = msh3trg->vertices_forming_triangles[i * 3];
		uint32_t n2 = msh3trg->vertices_forming_triangles[i*3+1];
		uint32_t n3 = msh3trg->vertices_forming_triangles[i*3+2];
		nb_graphics_move_to(g,
				    msh3trg->vertices[n1 * 2],
				    msh3trg->vertices[n1*2+1]);
		nb_graphics_line_to(g,
				    msh3trg->vertices[n2 * 2],
				    msh3trg->vertices[n2*2+1]);
		nb_graphics_line_to(g,
				    msh3trg->vertices[n3 * 2],
				    msh3trg->vertices[n3*2+1]);
		nb_graphics_close_path(g);

		nb_graphics_set_source_rgba(g, 25, 75, 255, 128);

		nb_graphics_fill_preserve(g);

		nb_graphics_set_source(g, NB_BLUE);

		nb_graphics_stroke(g);
	}

	/* Draw input segments */
	nb_graphics_set_source(g, NB_AQUAMARIN);

	nb_graphics_set_line_width(g, 1.0);
	for (uint32_t i = 0; i < msh3trg->N_input_segments; i++) {
		uint32_t N_vtx = msh3trg->N_vtx_x_inputsgm[i];
		if (0 < N_vtx) {
			uint32_t n1 = msh3trg->meshvtx_x_inputsgm[i][0];
			nb_graphics_move_to(g,
					    msh3trg->vertices[n1 * 2],
					    msh3trg->vertices[n1*2+1]);
			for (uint32_t j = 1; j < N_vtx; j++) {
				uint32_t n2 = msh3trg->meshvtx_x_inputsgm[i][j];
				nb_graphics_line_to(g,
						    msh3trg->vertices[n2 * 2],
						    msh3trg->vertices[n2*2+1]);
			}
		}
		nb_graphics_stroke(g);
	}
}

static void calculate_partition_centers
                     (const vcn_msh3trg_t *const restrict msh3trg,
		      const uint32_t *const restrict part,
		      uint32_t kpart, double **pcenter)
{
	for (uint32_t k = 0; k < kpart; k++)
		memset(pcenter[k], 0, 2 * sizeof(double));

	uint32_t* ksize = calloc(kpart, sizeof(*ksize));
	for (uint32_t i = 0; i < msh3trg->N_triangles; i++) {
		double xc = 0.0;
		double yc = 0.0;
		for (uint32_t j = 0; j < 3; j++) {
			uint32_t n = msh3trg->vertices_forming_triangles[i*3+j];
			xc += msh3trg->vertices[n * 2];
			yc += msh3trg->vertices[n*2+1];
		}
		pcenter[part[i]][0] += xc / 3.0;
		pcenter[part[i]][1] += yc / 3.0;
		ksize[part[i]] += 1;
	}

	for (uint32_t k = 0; k < kpart; k++) {
		pcenter[k][0] /= ksize[k];
		pcenter[k][1] /= ksize[k];
	}
}

static void draw_triangle_partition_edge(nb_graphics_context_t *g,
					 const vcn_msh3trg_t *const msh3trg,
					 const uint32_t *const part,
					 uint32_t kpart, double **pcenter,
					 double scale_partitions,
					 int trg_id, int edge_id)
{
	double vtx[2];
	uint32_t n1 = msh3trg->vertices_forming_triangles[trg_id*3+edge_id];
	uint32_t n2 = msh3trg->vertices_forming_triangles[trg_id*3 +
							  (edge_id+1)%3];

	memcpy(vtx, &(msh3trg->vertices[n1*2]), 2 * sizeof(double));
	if (scale_partitions > 0.0 && scale_partitions < 1.0)
		scale_vtx(vtx, pcenter[part[trg_id]], scale_partitions);

	nb_graphics_move_to(g, vtx[0], vtx[1]);

	memcpy(vtx, &(msh3trg->vertices[n2*2]), 2 * sizeof(double));
	if (scale_partitions > 0.0 && scale_partitions < 1.0)
		scale_vtx(vtx, pcenter[part[trg_id]], scale_partitions);

	nb_graphics_line_to(g, vtx[0], vtx[1]);
	nb_graphics_stroke(g);
}

static void draw_triangle_partition_border(nb_graphics_context_t *g,
					   const vcn_msh3trg_t *const msh3trg,
					   const uint32_t *const part,
					   uint32_t kpart, double **pcenter,
					   double scale_partitions,
					   uint32_t trg_id)
{
	uint32_t n1 = msh3trg->triangles_sharing_sides[trg_id * 3];
	uint32_t n2 = msh3trg->triangles_sharing_sides[trg_id*3+1];
	uint32_t n3 = msh3trg->triangles_sharing_sides[trg_id*3+2];

	if (part[trg_id] != part[n1])
		draw_triangle_partition_edge(g, msh3trg, part, kpart, pcenter,
					     scale_partitions, trg_id, 0);

	if (part[trg_id] != part[n2])
		draw_triangle_partition_edge(g, msh3trg, part, kpart, pcenter,
					     scale_partitions, trg_id, 1);

	if (part[trg_id] != part[n3])
		draw_triangle_partition_edge(g, msh3trg, part, kpart, pcenter,
					     scale_partitions, trg_id, 2);
}

static inline void scale_vtx(double vtx[2], double center[2], double zoom)
{
	vtx[0] = (vtx[0] - center[0]) * zoom + center[0];
	vtx[1] = (vtx[1] - center[1]) * zoom + center[1];
}

void vcn_msh3trg_partition_draw(const vcn_msh3trg_t *const msh3trg,
				const char* filename, int width, int height,
				uint32_t k_part, const uint32_t *const part,
				uint32_t k_to_draw, double scale_partitions)
{
	if (k_part < 2) {
		nb_graphics_export(filename, width, height,
				   draw_msh3trg, msh3trg);
	} else {
		part_data_t part_data;
		part_data.msh3trg = msh3trg;
		part_data.k_part = k_part;
		part_data.part = part;
		part_data.k_to_draw = k_to_draw;
		part_data.scale = scale_partitions;
		nb_graphics_export(filename, width, height,
				   draw_msh3trg_partition,
				   &part_data);
	}
}

static void draw_msh3trg_partition(nb_graphics_context_t *g,
				   int width, int height,
				   const void *const part_data)
{
	const part_data_t *const data = part_data;
	const vcn_msh3trg_t *const msh3trg = data->msh3trg;

	/* Compute cam.center and cam.zoom */
	double box[4];
	vcn_utils2D_get_enveloping_box_from_subset(msh3trg->N_input_vertices,
						   msh3trg->input_vertices,
						   msh3trg->vertices,
						   2 * sizeof(*(msh3trg->vertices)),
						   vcn_utils2D_get_x_from_darray,
						   vcn_utils2D_get_y_from_darray,
						   box);

	nb_graphics_enable_camera(g);
	nb_graphics_camera_t* cam = nb_graphics_get_camera(g);
	nb_graphics_cam_fit_box(cam, box, width, height);

	/* Calculate center of partitions */
	double** pcenter = NULL;
	if (data->scale > 0.0 && data->scale < 1.0) {
		pcenter = malloc(data->k_part * sizeof(*pcenter));
		for (uint32_t k = 0; k < data->k_part; k++)
			pcenter[k] = calloc(2, sizeof(*(pcenter[k])));
		calculate_partition_centers(msh3trg, data->part,
					    data->k_part, pcenter);
	}

	/* Draw triangles */
	for (uint32_t i = 0; i < msh3trg->N_triangles; i++) {
		if (data->part[i] != data->k_to_draw &&
		    data->k_to_draw < data->k_part)
			continue;
    
		uint32_t n1 = msh3trg->vertices_forming_triangles[i * 3];
		uint32_t n2 = msh3trg->vertices_forming_triangles[i*3+1];
		uint32_t n3 = msh3trg->vertices_forming_triangles[i*3+2];

		double vtx[2];
		memcpy(vtx, &(msh3trg->vertices[n1*2]), 2 * sizeof(double));
		if (data->scale > 0.0 && data->scale < 1.0)
			scale_vtx(vtx, pcenter[data->part[i]], data->scale);

		nb_graphics_move_to(g, vtx[0], vtx[1]);

		memcpy(vtx, &(msh3trg->vertices[n2*2]), 2 * sizeof(*vtx));
		if (data->scale > 0.0 && data->scale < 1.0)
			scale_vtx(vtx, pcenter[data->part[i]], data->scale);
		nb_graphics_line_to(g, vtx[0], vtx[1]);

		memcpy(vtx, &(msh3trg->vertices[n3*2]), 2 * sizeof(*vtx));
		if (data->scale > 0.0 && data->scale < 1.0)
			scale_vtx(vtx, pcenter[data->part[i]], data->scale);

		nb_graphics_line_to(g, vtx[0], vtx[1]);
		nb_graphics_close_path(g);

		nb_graphics_set_source(g, NB_LIGHT_GRAY);
		nb_graphics_fill_preserve(g);

		nb_graphics_set_source_rgb(g, 75, 128, 255);
		nb_graphics_set_line_width(g, 0.5);

		nb_graphics_stroke(g);

		/* Draw border */
		nb_graphics_set_source_rgb(g, 255, 0, 200);
		nb_graphics_set_line_width(g, 1.0);

		draw_triangle_partition_border(g, msh3trg, data->part,
					       data->k_part,
					       pcenter, data->scale, i);
	}

	if (data->scale > 0.0 && data->scale < 1.0) {
		for (uint32_t k = 0; k < data->k_part; k++)
			free(pcenter[k]);
		free(pcenter);
	}

	/* Draw input segments */
	if (data->k_to_draw >= data->k_part) {
		nb_graphics_set_source_rgb(g, 255, 0, 200);

		nb_graphics_set_line_width(g, 1.0);
		for (uint32_t i = 0; i < msh3trg->N_input_segments; i++) {
			uint32_t N_vtx = msh3trg->N_vtx_x_inputsgm[i];
			if (0 < N_vtx) {
				uint32_t n1 = msh3trg->meshvtx_x_inputsgm[i][0];
				nb_graphics_move_to(g,
						    msh3trg->vertices[n1 * 2],
						    msh3trg->vertices[n1*2+1]);
				for (uint32_t j = 0; j < N_vtx; j++) {
					uint32_t n2 = msh3trg->meshvtx_x_inputsgm[i][j];
					nb_graphics_line_to(g,
							    msh3trg->vertices[n2 * 2],
							    msh3trg->vertices[n2*2+1]);
				}
				nb_graphics_stroke(g);
			}
		}
	}
}
