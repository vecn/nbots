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
#include "vcn/pde_bot/finite_element/modules/exporter_cairo.h"

#include "../../../geometric_bot/mesh/mesh2D_structs.h"

static double _VCN_COLOR_BLACK[3] = {0.0, 0.0, 0.0};
static double _VCN_COLOR_BLUE[3] = {0.0, 0.0, 1.0};

typedef struct {
	double center[2];
	double zoom;
} camera_t;

static void set_center_and_zoom(camera_t *cam, double box[4],
				double width, double height);
static void draw_triangles(cairo_t *const cr,
			   const vcn_container_t *const ht_trg,
			   int width, int height,
			   const camera_t *const cam,
			   bool fill_triangles,
			   double rgba_bg[4], double rgb_fg[3]);
static void draw_input_vertices(cairo_t *const cr,
				uint32_t N_vertices,
				msh_vtx_t **vertices,
				int width, int height,
				const camera_t *const cam,
				bool include_numbering,
				double rgb_bg[3], double rgb_fg[3]);
static void draw_input_segments(cairo_t *const cr,
				uint32_t N_input_segments,
				msh_edge_t **segments,
				int width, int height,
				const camera_t *const cam,
				double rgb[3]);
static void mesh_draw_with_cairo(const vcn_mesh_t *const mesh,
				 cairo_t *cr, int width, int height,
				 const camera_t *const cam);
static double msh_vtx_get_x(const void *const vtx_ptr);
static double msh_vtx_get_y(const void *const vtx_ptr);
static void set_camera(camera_t *cam, const vcn_mesh_t *const mesh,
		       int width, int height);

static inline void set_center_and_zoom(camera_t *cam, double box[4],
				       double width, double height)
{
	cam->center[0] = (box[0] + box[2]) / 2.0;
	cam->center[1] = (box[1] + box[3]) / 2.0;
	cam->zoom = width / (box[2] - box[0]);
	if (cam->zoom > height / (box[3] - box[1]))
		cam->zoom = height / (box[3] - box[1]);
	cam->zoom *= 0.9;
}

static void draw_triangles(cairo_t *const cr,
			   const vcn_container_t *const ht_trg,
			   int width, int height,
			   const camera_t *const cam,
			   bool fill_triangles,
			   double rgba_bg[4], double rgb_fg[3])
{
	vcn_iterator_t* iter = vcn_iterator_create();
	vcn_iterator_set_container(iter, ht_trg);
	while (vcn_iterator_has_more(iter)) {
		const msh_trg_t* trg = vcn_iterator_get_next(iter);
		cairo_move_to(cr, 
			      cam->zoom * trg->v1->x[0] - 
			      cam->zoom * cam->center[0] + width / 2.0,
			      -cam->zoom * trg->v1->x[1] +
			      cam->zoom * cam->center[1] + height / 2.0);
		cairo_line_to(cr, 
			      cam->zoom * trg->v2->x[0] -
			      cam->zoom * cam->center[0] + width / 2.0,
			      -cam->zoom * trg->v2->x[1] +
			      cam->zoom * cam->center[1] + height / 2.0);
		cairo_line_to(cr, 
			      cam->zoom * trg->v3->x[0] -
			      cam->zoom * cam->center[0] + width / 2.0,
			      -cam->zoom * trg->v3->x[1] +
			      cam->zoom * cam->center[1] + height / 2.0);
		cairo_close_path(cr);
		if (fill_triangles) {
			cairo_set_source_rgba(cr, rgba_bg[0],
					      rgba_bg[1], rgba_bg[2],
					      rgba_bg[3]);
			cairo_fill_preserve(cr);
		}
		cairo_set_source_rgb(cr, rgb_fg[0], rgb_fg[1], rgb_fg[2]);
		cairo_stroke(cr);
	}
	vcn_iterator_destroy(iter);
}

static void draw_input_vertices(cairo_t *const restrict cr,
				uint32_t N_vertices,
				msh_vtx_t **vertices,
				int width, int height,
				const camera_t *const cam,
				bool include_numbering,
				double rgb_bg[3], double rgb_fg[3])
{
	for (uint32_t i = 0; i < N_vertices; i++) {
		cairo_set_source_rgb(cr, rgb_bg[0], rgb_bg[1], rgb_bg[2]);

		cairo_move_to(cr,
			      cam->zoom * vertices[i]->x[0] - 
			      cam->zoom * cam->center[0] + width / 2.0 + 3.0,
			      -cam->zoom * vertices[i]->x[1] + cam->zoom * cam->center[1] +
			      height / 2.0);
		cairo_arc(cr,
			  cam->zoom * vertices[i]->x[0] - 
			  cam->zoom * cam->center[0] + width / 2.0,
			  -cam->zoom * vertices[i]->x[1] + cam->zoom * cam->center[1] +
			  height / 2.0, 3.0, 0.0, 2.0 * VCN_MATH_PI);
		cairo_fill(cr);

		/* Draw labels */
		if (include_numbering) {
			cairo_set_source_rgb(cr, rgb_fg[0],
					     rgb_fg[1], rgb_fg[2]);

			/* Show id */
			char str_id[5];
			sprintf(str_id, "%i", i);
			cairo_select_font_face(cr, "Sans",
					       CAIRO_FONT_SLANT_NORMAL,
					       CAIRO_FONT_WEIGHT_NORMAL);
			cairo_set_font_size(cr, 9);
			cairo_move_to(cr, 
				      cam->zoom * vertices[i]->x[0] - 
				      cam->zoom * cam->center[0] + width/2.0 + 4.0,
				      -cam->zoom * vertices[i]->x[1] + 
				      cam->zoom * cam->center[1] + height/2.0 - 1.0);
			cairo_show_text(cr, str_id);
		}
	}
}

static void draw_input_segments(cairo_t *const cr,
				uint32_t N_input_segments,
				msh_edge_t **segments,
				int width, int height,
				const camera_t *const cam,
				double rgb[3])
{
	cairo_set_source_rgb(cr, rgb[0], rgb[1], rgb[2]);
	cairo_set_line_width(cr, 0.75);
	for (uint32_t i = 0; i < N_input_segments; i++) {
		if (NULL == segments[i])
			continue;

		msh_edge_t* sgm = segments[i];

		cairo_move_to(cr, 
			      cam->zoom * sgm->v1->x[0] - cam->zoom * cam->center[0] +
			      width/2.0,
			      -cam->zoom * sgm->v1->x[1] + cam->zoom * cam->center[1] +
			      height/2.0);
		msh_vtx_t* vtx = sgm->v2;
		while (NULL != sgm) {
			cairo_line_to(cr, 
				      cam->zoom * vtx->x[0] - cam->zoom * cam->center[0] +
				      width/2.0,
				      -cam->zoom * vtx->x[1] + cam->zoom * cam->center[1] +
				      height/2.0);
			msh_edge_t* prev_sgm = sgm;
			sgm = medge_subsgm_next(sgm);
			if (NULL != sgm) {
				if (sgm->v1 == vtx) {
					vtx = sgm->v2;
				} else {
					if ((sgm->v1 == prev_sgm->v2) ||
					    (sgm->v1 == prev_sgm->v1))
						vtx = sgm->v2;
					else
						vtx = sgm->v1;
				}
			} else {
				break;
			}
		}
		cairo_stroke(cr);
	}
}

static void mesh_draw_with_cairo(const vcn_mesh_t *const mesh,
				 cairo_t* cr, int width, int height,
				 const camera_t *const cam)
{
	/* Draw background */
	cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
	cairo_paint(cr);

	/* Set back-ground and fore-ground colors */
	double rgba_bg[4] = {0.1, 0.1, 1.0, 0.5};

	/* Draw triangles */
	cairo_set_line_width(cr, 0.5);
	draw_triangles(cr, mesh->ht_trg, width, height,
		       cam, true, rgba_bg, 
		       _VCN_COLOR_BLUE);

	/* Draw input segments */
	if (true) {
		double rgb_segments[3] = {1.0, 0.0, 0.8};
		draw_input_segments(cr, mesh->N_input_sgm, mesh->input_sgm,
				    width, height, cam, rgb_segments);
	}
	/* Draw input vertices */
	if (true) {
		draw_input_vertices(cr, mesh->N_input_vtx, mesh->input_vtx,
				    width, height, cam, false,
				    _VCN_COLOR_BLUE, _VCN_COLOR_BLACK);
	}
}

static inline double msh_vtx_get_x(const void *const vtx_ptr)
{
	const msh_vtx_t *const *vtx = vtx_ptr;
	return (*vtx)->x[0];
}

static inline double msh_vtx_get_y(const void *const vtx_ptr)
{
	const msh_vtx_t *const *vtx = vtx_ptr;
	return (*vtx)->x[1];
}

void nb_fem_save_png(const vcn_mesh_t *const restrict mesh,
		     const double *results,
		     const char* filename,
		     int width, int height)
{
	cairo_surface_t* surface =
		cairo_image_surface_create(CAIRO_FORMAT_ARGB32, 
					   width, height);
	cairo_t* cr = cairo_create(surface);

	camera_t cam;
	set_camera(&cam, mesh, width, height);

	mesh_draw_with_cairo(mesh, cr, width, height, &cam);

	cairo_surface_write_to_png(surface, filename);

	cairo_destroy(cr);
	cairo_surface_destroy(surface);
}

static void set_camera(camera_t *cam, const vcn_mesh_t *const mesh,
		       int width, int height)
{
	double box[4];
	vcn_utils2D_get_enveloping_box(mesh->N_input_vtx, mesh->input_vtx,
				       sizeof(*(mesh->input_vtx)),
				       msh_vtx_get_x, msh_vtx_get_y,
				       box);
	set_center_and_zoom(cam, box, width, height);
}

void nb_fem_save_eps(const vcn_mesh_t *const mesh,
		     const double *results,
		     const char* filename,
		     int width, int height)
{
	/* Create drawable surface and cairo context */
	cairo_surface_t* surface = 
		cairo_ps_surface_create (filename, width, height);
	cairo_ps_surface_set_eps(surface, 1 /* TRUE from cairo_bool_t*/);

	cairo_t* cr = cairo_create(surface);

	/* Initialize Post script */
	cairo_ps_surface_dsc_begin_page_setup(surface);

	if(height < width)
		cairo_ps_surface_dsc_comment(surface,
					     "%%PageOrientation: Portrait");
	else
		cairo_ps_surface_dsc_comment(surface,
					     "%%PageOrientation: Landscape");

	camera_t cam;
	set_camera(&cam, mesh, width, height);

	mesh_draw_with_cairo(mesh, cr, width, height, &cam);

	cairo_surface_show_page(surface);/* Show post-script */

	cairo_destroy(cr);
	cairo_surface_finish(surface);
	cairo_surface_destroy(surface);
}
