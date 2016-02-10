#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <cairo.h>
#include <cairo-ps.h>

#include "vcn/math_bot.h"
#include "vcn/visual_cat.h"
#include "vcn/container_bot.h"
#include "vcn/geometric_bot.h"
#include "vcn/pde_bot/finite_element/modules/exporter_cairo.h"

typedef struct {
	double center[2];
	double zoom;
} camera_t;

static void set_center_and_zoom(camera_t *cam, double box[4],
				double width, double height);
static void set_camera_vtx(double v_dest[2], const double v_src[2],
			   camera_t *cam, double width, double height);
static void set_pattern_patch(cairo_pattern_t *pat, const double v1[2],
			      const double v2[2], const double v3[2]);
static void set_pattern_color(cairo_pattern_t *pat,
			      const vcn_palette_t *palette,
			      const double *results,
			      uint32_t min_id, uint32_t max_id,
			      uint32_t id1, uint32_t id2, uint32_t id3);

void nb_fem_save_png(const vcn_msh3trg_t *const msh3trg,
		     const double *results,
		     const char* filename,
		     int width, int height)
{

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
	set_center_and_zoom(&cam, box, width, height);

	/* Create drawable surface and cairo context */
	cairo_surface_t* surface =
		cairo_image_surface_create(CAIRO_FORMAT_ARGB32, width, height);
	cairo_t *const restrict cr = cairo_create(surface);

	/* Draw background */
	cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
	cairo_paint(cr);

	vcn_palette_t *palette = vcn_palette_create_preset(VCN_PALETTE_RAINBOW);

	uint32_t min_id;
	uint32_t max_id;
	vcn_array_get_min_max_ids(results, msh3trg->N_vertices, sizeof(*results),
				  vcn_compare_double, &min_id, &max_id);

	/* Draw triangles */
	cairo_set_line_width(cr, 0.5);
	for (uint32_t i = 0; i < msh3trg->N_triangles; i++) {
		uint32_t n1 = msh3trg->vertices_forming_triangles[i * 3];
		uint32_t n2 = msh3trg->vertices_forming_triangles[i*3+1];
		uint32_t n3 = msh3trg->vertices_forming_triangles[i*3+2];
		double v1[2];
		set_camera_vtx(v1, &(msh3trg->vertices[n1*2]), &cam,
			       width, height);
		double v2[2];
		set_camera_vtx(v2, &(msh3trg->vertices[n2*2]), &cam,
			       width, height);
		double v3[2];
		set_camera_vtx(v3, &(msh3trg->vertices[n3*2]), &cam,
			       width, height);
		cairo_move_to(cr, v1[0], v1[1]);
		cairo_line_to(cr, v2[0], v2[1]);
		cairo_line_to(cr, v3[0], v3[1]);
		cairo_close_path(cr);

		cairo_pattern_t *pat = cairo_pattern_create_mesh();
		set_pattern_patch(pat, v1, v2, v3);
		set_pattern_color(pat, palette, results,
				  min_id, max_id, n1, n2, n3);
		cairo_mesh_pattern_end_patch(pat);
		cairo_set_source(cr, pat);
		cairo_fill_preserve(cr);
		cairo_pattern_destroy(pat);

		cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
		cairo_stroke(cr);
	}

	/* Draw input segments */
	cairo_set_source_rgb(cr, 0.9, 0.0, 0.8);

	cairo_set_line_width(cr, 1.0);
	for (uint32_t i = 0; i < msh3trg->N_input_segments; i++) {
		if (msh3trg->N_subsgm_x_inputsgm[i] == 0)
			continue;

		uint32_t n1 = msh3trg->meshvtx_x_inputsgm[i][0];
		cairo_move_to(cr, 
			      cam.zoom * msh3trg->vertices[n1 * 2] -
			      cam.zoom * cam.center[0] + width / 2.0,
			      -cam.zoom * msh3trg->vertices[n1*2+1] +
			      cam.zoom * cam.center[1] + height / 2.0);
		for (uint32_t j = 0; j < msh3trg->N_subsgm_x_inputsgm[i]; j++) {
			uint32_t n2 = msh3trg->meshvtx_x_inputsgm[i][j+1];
			cairo_line_to(cr, 
				      cam.zoom * msh3trg->vertices[n2 * 2] -
				      cam.zoom * cam.center[0] + width/2.0,
				      -cam.zoom * msh3trg->vertices[n2*2+1] +
				      cam.zoom * cam.center[1] + height/2.0);
		}
		cairo_stroke(cr);
	}

	/* Write file */
	cairo_surface_write_to_png(surface, filename);
	/* Free cairo structures */
	cairo_destroy(cr);
	cairo_surface_destroy(surface);
}

static void set_center_and_zoom(camera_t *cam, double box[4],
				double width, double height)
{
	cam->center[0] = (box[0] + box[2]) / 2.0;
	cam->center[1] = (box[1] + box[3]) / 2.0;
	cam->zoom = width / (box[2] - box[0]);
	if (cam->zoom > height / (box[3] - box[1]))
		cam->zoom = height / (box[3] - box[1]);
	cam->zoom *= 0.9;
}

static void set_camera_vtx(double v_dest[2], const double v_src[2],
			   camera_t *cam, double width, double height)
{
	v_dest[0] = cam->zoom * (v_src[0] - cam->center[0]) + width/2.0;
	v_dest[1] = -cam->zoom * (v_src[1] - cam->center[1]) + height/2.0;
}

static void set_pattern_patch(cairo_pattern_t *pat, const double v1[2],
			      const double v2[2], const double v3[2])
{
	cairo_mesh_pattern_begin_patch(pat);
	cairo_mesh_pattern_move_to(pat, v1[0], v1[1]);
	cairo_mesh_pattern_line_to(pat, v2[0], v2[1]);
	cairo_mesh_pattern_line_to(pat, v3[0], v3[1]);
}

static void set_pattern_color(cairo_pattern_t *pat,
			      const vcn_palette_t *palette,
			      const double *results,
			      uint32_t min_id, uint32_t max_id,
			      uint32_t id1, uint32_t id2, uint32_t id3)
{
	double results_range = results[max_id] - results[min_id];
	uint8_t rgb[3];
	double val = (results[id1] - results[min_id]) / results_range;
	vcn_palette_get_colour(palette, val, rgb);
	cairo_mesh_pattern_set_corner_color_rgb(pat, 0, 
						rgb[0]/255.0f, 
						rgb[1]/255.0f, 
						rgb[2]/255.0f);
	val = (results[id2] - results[min_id]) / results_range;
	vcn_palette_get_colour(palette, val, rgb);
	cairo_mesh_pattern_set_corner_color_rgb(pat, 1, 
						rgb[0]/255.0f,
						rgb[1]/255.0f, 
						rgb[2]/255.0f);
	val = (results[id3] - results[min_id]) / results_range;
	vcn_palette_get_colour(palette, val, rgb);
	cairo_mesh_pattern_set_corner_color_rgb(pat, 2, 
						rgb[0]/255.0f, 
						rgb[1]/255.0f,
						rgb[2]/255.0f);
}
