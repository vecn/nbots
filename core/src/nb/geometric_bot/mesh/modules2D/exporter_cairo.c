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

#include "../mesh2D_structs.h"
#include "exporter_cairo/drawing_utils.h"

static void calculate_partition_centers(const vcn_msh3trg_t *const msh3trg,
					const uint32_t *const part,
					uint32_t kpart, double **pcenter);
static void draw_triangle_partition_edge(cairo_t *cr,
					 const vcn_msh3trg_t *const msh3trg,
					 const uint32_t *const part,
					 uint32_t kpart, double **pcenter,
					 double scale_partitions,
					 int trg_id, int edge_id,
					 const camera_t *const cam,
					 double width, double height);
static void draw_triangle_partition_border(cairo_t *cr,
					   const vcn_msh3trg_t *const msh3trg,
					   const uint32_t *const part,
					   uint32_t kpart, double **pcenter,
					   double scale_partitions,
					   uint32_t trg_id,
					   const camera_t *const cam,
					   double width, double height);
static void scale_vtx(double vtx[2], double center[2], double zoom);
static void draw_triangles(cairo_t *const cr,
			   const nb_container_t *const ht_trg,
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
static void draw_vertices_halos(cairo_t *cr, 
				uint32_t N, msh_vtx_t **vertices,
				int width, int height,
				const camera_t *const cam);
static void set_camera(camera_t *cam, const vcn_mesh_t *const mesh,
		       int width, int height);

static void draw_triangles(cairo_t *const cr,
			   const nb_container_t *const ht_trg,
			   int width, int height,
			   const camera_t *const cam,
			   bool fill_triangles,
			   double rgba_bg[4], double rgb_fg[3])
{
	nb_iterator_t* iter = nb_iterator_create();
	nb_iterator_set_container(iter, ht_trg);
	while (nb_iterator_has_more(iter)) {
		const msh_trg_t* trg = nb_iterator_get_next(iter);
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
	nb_iterator_destroy(iter);
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
			  height / 2.0, 3.0, 0.0, 2.0 * NB_MATH_PI);
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
		       _NB_RGB_BLUE);

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
				    _NB_RGB_BLUE, _NB_RGB_BLACK);
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

void vcn_mesh_save_png(const vcn_mesh_t *const restrict mesh,
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

void vcn_dewall_save_png(const vcn_mesh_t *const restrict mesh,
			 const char* filename, int width, int height,
			 uint8_t axe, double alpha, uint32_t N,
			 void *vtx_array)
{
	msh_vtx_t **vertices = vtx_array;

	/* Create drawable surface and cairo context */
	cairo_surface_t* surface =
		cairo_image_surface_create(CAIRO_FORMAT_ARGB32, 
					   width, height);
	cairo_t* cr = cairo_create(surface);

	camera_t cam;
	set_camera(&cam, mesh, width, height);

	mesh_draw_with_cairo(mesh, cr, width, height, &cam);

	/* Draw Wall */
	cairo_set_source_rgb(cr, 1.0, 0.0, 0.0);
	if (0 == axe) {
		double x_alpha = cam.zoom * alpha -
			cam.zoom * cam.center[0] + width / 2.0;
		cairo_move_to(cr, x_alpha, 0);
		cairo_line_to(cr, x_alpha, height);
	} else {
		double y_alpha = -cam.zoom * alpha +
			cam.zoom * cam.center[1] + height / 2.0;
		cairo_move_to(cr, 0, y_alpha);
		cairo_line_to(cr, width, y_alpha);
	}
	cairo_stroke(cr);

	/* Draw vertices of partition */
	cairo_set_source_rgb(cr, 0.0, 1.0, 0.0);
	cairo_set_line_width(cr, 1.0);
	draw_vertices_halos(cr, N, vertices, width, height, &cam);

	/* Write file */
	cairo_surface_write_to_png(surface, filename);

	/* Free cairo structures */
	cairo_destroy(cr);
	cairo_surface_destroy(surface);
}

static void draw_vertices_halos(cairo_t *cr,
				uint32_t N, msh_vtx_t **vertices,
				int width, int height,
				const camera_t *const cam)
{
	for (uint32_t i = 0; i < N; i++) {
		cairo_move_to(cr,
			      cam->zoom * vertices[i]->x[0] - 
			      cam->zoom * cam->center[0] + width / 2.0 + 3.0,
			      -cam->zoom * vertices[i]->x[1] +
			      cam->zoom * cam->center[1] + height / 2.0);
		cairo_arc(cr,
			  cam->zoom * vertices[i]->x[0] - 
			  cam->zoom * cam->center[0] + width / 2.0,
			  -cam->zoom * vertices[i]->x[1] + 
			  cam->zoom * cam->center[1] +
			  height / 2.0, 5, 0.0, 2.0 * NB_MATH_PI);
		cairo_stroke(cr);
	}
}

static void set_camera(camera_t *cam, const vcn_mesh_t *const mesh,
		       int width, int height)
{
	double box[4];
	vcn_utils2D_get_enveloping_box(mesh->N_input_vtx, mesh->input_vtx,
				       sizeof(*(mesh->input_vtx)),
				       msh_vtx_get_x, msh_vtx_get_y,
				       box);
	nb_drawing_utils_set_center_and_zoom(cam, box, width, height);
}

void vcn_mesh_save_eps(const vcn_mesh_t *const mesh,
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

void vcn_msh3trg_save_png(const vcn_msh3trg_t *const restrict msh3trg,
			  const char* filename, int width, int height,
			  double rgba_bg[4], double rgba_fg[4],
			  double rgb_edge[3], double rgb_sgm[3],
			  double edge_width, double sgm_width)
{
	/* Compute cam->center and cam->zoom */
	double box[4];
	vcn_utils2D_get_enveloping_box_from_subset(msh3trg->N_input_vertices,
						   msh3trg->input_vertices,
						   msh3trg->vertices,
						   2 * sizeof(*(msh3trg->vertices)),
						   vcn_utils2D_get_x_from_darray,
						   vcn_utils2D_get_y_from_darray,
						   box);

	camera_t cam;
	nb_drawing_utils_set_center_and_zoom(&cam, box, width, height);

	/* Create drawable surface and cairo context */
	cairo_surface_t* surface =
		cairo_image_surface_create(CAIRO_FORMAT_ARGB32, width, height);
	cairo_t *const restrict cr = cairo_create(surface);

	/* Draw background */
	if (NULL == rgba_bg)
		cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
	else 
		cairo_set_source_rgba(cr, rgba_bg[0], rgba_bg[1],
				      rgba_bg[2], rgba_bg[3]);
	cairo_paint(cr);

	/* Draw triangles */
	cairo_set_line_width(cr, edge_width);
	for (uint32_t i = 0; i < msh3trg->N_triangles; i++) {
		uint32_t n1 = msh3trg->vertices_forming_triangles[i * 3];
		uint32_t n2 = msh3trg->vertices_forming_triangles[i*3+1];
		uint32_t n3 = msh3trg->vertices_forming_triangles[i*3+2];
		cairo_move_to(cr, 
			      cam.zoom * msh3trg->vertices[n1 * 2] -
			      cam.zoom * cam.center[0] + width/2.0,
			      -cam.zoom * msh3trg->vertices[n1*2+1] +
			      cam.zoom * cam.center[1] + height/2.0);
		cairo_line_to(cr, 
			      cam.zoom * msh3trg->vertices[n2 * 2] -
			      cam.zoom * cam.center[0] + width/2.0,
			      -cam.zoom * msh3trg->vertices[n2*2+1] +
			      cam.zoom * cam.center[1] + height/2.0);
		cairo_line_to(cr, 
			      cam.zoom * msh3trg->vertices[n3 * 2] -
			      cam.zoom * cam.center[0] + width/2.0,
			      -cam.zoom * msh3trg->vertices[n3*2+1] +
			      cam.zoom * cam.center[1] + height/2.0);
		cairo_close_path(cr);

		if (NULL == rgba_fg)
			cairo_set_source_rgba(cr, 0.3, 0.5, 1.0, 0.25);
		else
			cairo_set_source_rgba(cr, rgba_fg[0], rgba_fg[1],
					      rgba_fg[2], rgba_fg[3]);
		cairo_fill_preserve(cr);

		if (NULL == rgb_edge)
			cairo_set_source_rgb(cr, 0.0, 0.3, 1.0);
		else
			cairo_set_source_rgb(cr, rgb_edge[0],
					     rgb_edge[1], rgb_edge[2]);
		cairo_stroke(cr);
	}

	/* Draw input segments */
	if (NULL == rgb_sgm)
		cairo_set_source_rgb(cr, 0.9, 0.0, 0.8);
	else
		cairo_set_source_rgb(cr, rgb_sgm[0], rgb_sgm[1], rgb_sgm[2]);

	cairo_set_line_width(cr, sgm_width);
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

static void draw_triangle_partition_edge(cairo_t *cr,
					 const vcn_msh3trg_t *const msh3trg,
					 const uint32_t *const part,
					 uint32_t kpart, double **pcenter,
					 double scale_partitions,
					 int trg_id, int edge_id,
					 const camera_t *const cam,
					 double width, double height)
{
	double vtx[2];
	uint32_t n1 = msh3trg->vertices_forming_triangles[trg_id*3+edge_id];
	uint32_t n2 = msh3trg->vertices_forming_triangles[trg_id*3+(edge_id+1)%3];

	memcpy(vtx, &(msh3trg->vertices[n1*2]), 2 * sizeof(double));
	if (scale_partitions > 0.0 && scale_partitions < 1.0)
		scale_vtx(vtx, pcenter[part[trg_id]], scale_partitions);
	cairo_move_to(cr, cam->zoom * vtx[0] -
		      cam->zoom * cam->center[0] + width/2.0,
		      -cam->zoom * vtx[1] + cam->zoom * 
		      cam->center[1] + height/2.0);

	memcpy(vtx, &(msh3trg->vertices[n2*2]), 2 * sizeof(double));
	if (scale_partitions > 0.0 && scale_partitions < 1.0)
		scale_vtx(vtx, pcenter[part[trg_id]], scale_partitions);
	cairo_line_to(cr, cam->zoom * vtx[0] -
		      cam->zoom * cam->center[0] + width/2.0,
		      -cam->zoom * vtx[1] + cam->zoom * 
		      cam->center[1] + height/2.0);
	cairo_stroke(cr);
}

static void draw_triangle_partition_border(cairo_t *cr,
					   const vcn_msh3trg_t *const msh3trg,
					   const uint32_t *const part,
					   uint32_t kpart, double **pcenter,
					   double scale_partitions,
					   uint32_t trg_id,
					   const camera_t *const cam,
					   double width, double height)
{
	uint32_t n1 = msh3trg->triangles_sharing_sides[trg_id * 3];
	uint32_t n2 = msh3trg->triangles_sharing_sides[trg_id*3+1];
	uint32_t n3 = msh3trg->triangles_sharing_sides[trg_id*3+2];

	if (part[trg_id] != part[n1])
		draw_triangle_partition_edge(cr, msh3trg, part, kpart, pcenter,
					     scale_partitions, trg_id, 0,
					     cam, width, height);

	if (part[trg_id] != part[n2])
		draw_triangle_partition_edge(cr, msh3trg, part, kpart, pcenter,
					     scale_partitions, trg_id, 1,
					     cam, width, height);

	if (part[trg_id] != part[n3])
		draw_triangle_partition_edge(cr, msh3trg, part, kpart, pcenter,
					     scale_partitions, trg_id, 2,
					     cam, width, height);
}

static inline void scale_vtx(double vtx[2], double center[2], double zoom)
{
	vtx[0] = (vtx[0] - center[0]) * zoom + center[0];
	vtx[1] = (vtx[1] - center[1]) * zoom + center[1];
}

void vcn_msh3trg_partition_save_png(const vcn_msh3trg_t *const msh3trg,
				    const char* filename, int width, int height,
				    uint32_t k_part, const uint32_t *const part,
				    uint32_t k_to_draw, double scale_partitions,
				    double rgba_bg[4], double alpha_fg,
				    double rgb_edge[3], double rgb_sgm[3],
				    double edge_width, double sgm_width,
				    bool use_colors)
{
	if (k_part < 2)
		vcn_msh3trg_save_png(msh3trg, filename, width, height,
				     rgba_bg, NULL, rgb_edge, rgb_sgm,
				     edge_width, sgm_width);

	/* Define colors */
	uint32_t N_part_colors = 8;
	double part_color[8][3] = {{1.0, 0.0, 0.0},
				   {1.0, 1.0, 0.0},
				   {0.0, 0.0, 1.0},
				   {0.0, 1.0, 0.0},
				   {0.5, 0.0, 1.0},
				   {1.0, 0.5, 0.0},
				   {0.0, 1.0, 1.0},
				   {1.0, 0.0, 1.0}};

	/* Compute cam.center and cam.zoom */
	double box[4];
	vcn_utils2D_get_enveloping_box_from_subset(msh3trg->N_input_vertices,
						   msh3trg->input_vertices,
						   msh3trg->vertices,
						   2 * sizeof(*(msh3trg->vertices)),
						   vcn_utils2D_get_x_from_darray,
						   vcn_utils2D_get_y_from_darray,
						   box);

	camera_t cam;
	nb_drawing_utils_set_center_and_zoom(&cam, box, width, height);

	/* Calculate cam.center of partitions */
	double** pcenter = NULL;
	if (scale_partitions > 0.0 && scale_partitions < 1.0) {
		pcenter = malloc(k_part * sizeof(*pcenter));
		for (uint32_t k = 0; k < k_part; k++)
			pcenter[k] = calloc(2, sizeof(*(pcenter[k])));
		calculate_partition_centers(msh3trg, part, k_part, pcenter);
	}

	/* Create drawable surface and cairo context */
	cairo_surface_t* surface =
		cairo_image_surface_create(CAIRO_FORMAT_ARGB32, width, height);
	cairo_t *const restrict cr = cairo_create(surface);

	/* Draw background */
	if (NULL == rgba_bg)
		cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
	else 
		cairo_set_source_rgba(cr, rgba_bg[0], rgba_bg[1], rgba_bg[2], rgba_bg[3]);
	cairo_paint(cr);

	/* Draw triangles */
	for (uint32_t i = 0; i < msh3trg->N_triangles; i++) {
		if (part[i] != k_to_draw && k_to_draw < k_part)
			continue;
    
		uint32_t n1 = msh3trg->vertices_forming_triangles[i * 3];
		uint32_t n2 = msh3trg->vertices_forming_triangles[i*3+1];
		uint32_t n3 = msh3trg->vertices_forming_triangles[i*3+2];

		double vtx[2];
		memcpy(vtx, &(msh3trg->vertices[n1*2]), 2 * sizeof(double));
		if (scale_partitions > 0.0 && scale_partitions < 1.0)
			scale_vtx(vtx, pcenter[part[i]], scale_partitions);

		cairo_move_to(cr, cam.zoom * vtx[0] - cam.zoom * cam.center[0] + width/2.0,
			      -cam.zoom * vtx[1] + cam.zoom * cam.center[1] + height/2.0);

		memcpy(vtx, &(msh3trg->vertices[n2*2]), 2 * sizeof(*vtx));
		if (scale_partitions > 0.0 && scale_partitions < 1.0)
			scale_vtx(vtx, pcenter[part[i]], scale_partitions);
		cairo_line_to(cr, cam.zoom * vtx[0] - cam.zoom * cam.center[0] + width/2.0,
			      -cam.zoom * vtx[1] + cam.zoom * cam.center[1] + height/2.0);

		memcpy(vtx, &(msh3trg->vertices[n3*2]), 2 * sizeof(*vtx));
		if (scale_partitions > 0.0 && scale_partitions < 1.0)
			scale_vtx(vtx, pcenter[part[i]], scale_partitions);
		cairo_line_to(cr, cam.zoom * vtx[0] - cam.zoom * cam.center[0] + width/2.0,
			      -cam.zoom * vtx[1] + cam.zoom * cam.center[1] + height/2.0);
		cairo_close_path(cr);

		uint32_t cid = part[i] % N_part_colors;
		if (use_colors)
			cairo_set_source_rgba(cr, part_color[cid][0],
					      part_color[cid][1],
					      part_color[cid][2], alpha_fg);
		else
			cairo_set_source_rgba(cr, 1.0, 1.0, 1.0, alpha_fg);
		cairo_fill_preserve(cr);

		if (NULL == rgb_edge)
			cairo_set_source_rgb(cr, 0.3, 0.5, 1.0);
		else
			cairo_set_source_rgb(cr, rgb_edge[0], rgb_edge[1], rgb_edge[2]);
		cairo_set_line_width(cr, edge_width);

		cairo_stroke(cr);

		/* Draw border */
		if (NULL == rgb_sgm)
			cairo_set_source_rgb(cr, 1.0, 0.0, 0.8);
		else
			cairo_set_source_rgb(cr, rgb_sgm[0], rgb_sgm[1], rgb_sgm[2]);
		cairo_set_line_width(cr, sgm_width);

		draw_triangle_partition_border(cr, msh3trg, part, k_part,
					       pcenter, scale_partitions, i,
					       &cam, width, height);
	}

	if (scale_partitions > 0.0 && scale_partitions < 1.0) {
		for (uint32_t k = 0; k < k_part; k++)
			free(pcenter[k]);
		free(pcenter);
	}

	/* Draw input segments */
	if (k_to_draw >= k_part) {
		if (NULL == rgb_sgm)
			cairo_set_source_rgb(cr, 1.0, 0.0, 0.8);
		else
			cairo_set_source_rgb(cr, rgb_sgm[0],
					     rgb_sgm[1], rgb_sgm[2]);

		cairo_set_line_width(cr, sgm_width);
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
	}

	/* Write file */
	cairo_surface_write_to_png(surface, filename);
	/* Free cairo structures */
	cairo_destroy(cr);
	cairo_surface_destroy(surface);
}

void vcn_mshpack_save_png(const vcn_mshpack_t *const restrict mshpack,
			  const char* filename,
			  int width, int height)
{
	if (0 == mshpack->N_spheres)
		return;

	/* Compute cam.center and cam.zoom */
	double xmin = mshpack->centers[0] - mshpack->radii[0];
	double xmax = mshpack->centers[0] + mshpack->radii[0];
	double ymin = mshpack->centers[1] - mshpack->radii[0];
	double ymax = mshpack->centers[1] + mshpack->radii[0];
	for (uint32_t i = 1; i < mshpack->N_spheres; i++) {
		if (xmin > mshpack->centers[i * 2] - mshpack->radii[i])
			xmin = mshpack->centers[i * 2] - mshpack->radii[i];
		else if (xmax < mshpack->centers[i * 2] + mshpack->radii[i]) 
			xmax = mshpack->centers[i * 2] + mshpack->radii[i];
		if (ymin > mshpack->centers[i*2+1] - mshpack->radii[i]) 
			ymin = mshpack->centers[i*2+1] - mshpack->radii[i];
		else if (ymax < mshpack->centers[i*2+1] + mshpack->radii[i])
			ymax = mshpack->centers[i*2+1] + mshpack->radii[i];
	}
	camera_t cam;
	cam.center[0] = (xmin+xmax)/2.0;
	cam.center[1] = (ymin+ymax)/2.0;
	cam.zoom = width/(xmax-xmin);
	if (cam.zoom > height/(ymax-ymin))
		cam.zoom = height/(ymax-ymin);
	cam.zoom *= 0.9;

	/* Create drawable surface and cairo context */
	cairo_surface_t* surface =
		cairo_image_surface_create(CAIRO_FORMAT_ARGB32, 
					   width, height);
	cairo_t *const restrict cr = cairo_create(surface);

	/* Draw background */
	cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
	cairo_paint(cr);

	/* Draw spheres */
	for (uint32_t i = 0; i < mshpack->N_spheres; i++) {
		cairo_arc(cr, 
			  cam.zoom * mshpack->centers[i * 2] -
			  cam.zoom * cam.center[0] + width/2.0,
			  -cam.zoom * mshpack->centers[i*2+1] +
			  cam.zoom * cam.center[1] + height/2.0,
			  mshpack->radii[i] * cam.zoom, 0, 2 * NB_MATH_PI);
		if (true) {
			cairo_set_source_rgba(cr, 0.3, 0.5, 1.0, 0.35);
			cairo_fill_preserve(cr);
		}
		cairo_set_line_width(cr, 0.5);
		cairo_set_source_rgb(cr, 0.3, 0.5, 1.0);
		cairo_stroke(cr);
	}

	/* Write file */
	cairo_surface_write_to_png(surface, filename);
	/* Free cairo structures */
	cairo_destroy(cr);
	cairo_surface_destroy(surface);
}
