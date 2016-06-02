#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>

#include "nb/math_bot.h"
#include "nb/container_bot.h"
#include "nb/geometric_bot.h"
#include "nb/graphics_bot.h"

#include "../../mesh2D_structs.h"

typedef struct {
	const nb_mesh_t *mesh;
	uint8_t axe;
	double alpha;
	uint32_t N;
	void *vtx_array;	
} dewall_data_t;

static void draw_triangles(void *const draw_ptr,
			   const nb_container_t *const ht_trg,
			   int width, int height,
			   const camera_t *const cam);
static void draw_input_segments(void *const draw_ptr,
				uint32_t N_input_segments,
				msh_edge_t **segments,
				int width, int height,
				const camera_t *const cam);
static void draw_mesh(void *draw_ptr, int width, int height,
		      const void *const mesh_ptr);
static double msh_vtx_get_x(const void *const vtx_ptr);
static double msh_vtx_get_y(const void *const vtx_ptr);

static void draw_dewall(void *draw_ptr, int width, int height,
			const void *const mesh_ptr);
static void draw_vertices_halos(void *draw_ptr, 
				uint32_t N, msh_vtx_t **vertices,
				int width, int height,
				const camera_t *const cam);


static void draw_triangles(void *const draw_ptr,
			   const nb_container_t *const ht_trg,
			   int width, int height,
			   const camera_t *const cam)
{
	nb_iterator_t* iter = nb_iterator_create();
	nb_iterator_set_container(iter, ht_trg);
	while (nb_iterator_has_more(iter)) {
		const msh_trg_t* trg = nb_iterator_get_next(iter);
		nb_drawing_move_to(draw_ptr, cam, trg->v1->x[0], trg->v1->x[1]);
		nb_drawing_line_to(draw_ptr, cam, trg->v2->x[0], trg->v2->x[1]);
		nb_drawing_line_to(draw_ptr, cam, trg->v3->x[0], trg->v3->x[1]);
		nb_drawing_close_path(draw_ptr);
		nb_drawing_set_source_rgba(draw_ptr, 0.1,
					   0.3, 1.0, 0.5);
		nb_drawing_fill_preserve(draw_ptr);
		nb_drawing_set_source_rgb(draw_ptr, 0, 0, 1);
		nb_drawing_stroke(draw_ptr);
	}
	nb_iterator_destroy(iter);
}

static void draw_input_segments(void *const draw_ptr,
				uint32_t N_input_segments,
				msh_edge_t **segments,
				int width, int height,
				const camera_t *const cam)
{
  nb_drawing_set_source_rgb(draw_ptr, 0.9, 0.2, 0.4);
	nb_drawing_set_line_width(draw_ptr, 0.75);
	for (uint32_t i = 0; i < N_input_segments; i++) {
		if (NULL == segments[i])
			continue;

		msh_edge_t* sgm = segments[i];

		nb_drawing_move_to(draw_ptr, cam, sgm->v1->x[0],
				   sgm->v1->x[1]);
		msh_vtx_t* vtx = sgm->v2;
		while (NULL != sgm) {
			nb_drawing_line_to(draw_ptr, cam,
					   vtx->x[0], vtx->x[1]);
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
		nb_drawing_stroke(draw_ptr);
	}
}

static void draw_mesh(void *draw_ptr, int width, int height,
		      const void *const mesh_ptr)
{
	const nb_mesh_t *const mesh = mesh_ptr;

	double box[4];
	vcn_utils2D_get_enveloping_box(mesh->N_input_vtx,
				       mesh->input_vtx,
				       sizeof(*(mesh->input_vtx)),
				       msh_vtx_get_x,
				       msh_vtx_get_y,
				       box);

	camera_t cam;
	nb_drawing_utils_set_center_and_zoom(&cam, box, width, height);

	/* Draw triangles */
	nb_drawing_set_line_width(draw_ptr, 0.5);
	draw_triangles(draw_ptr, mesh->ht_trg, width, height, &cam);

	/* Draw input segments */
	draw_input_segments(draw_ptr, mesh->N_input_sgm, mesh->input_sgm,
			    width, height, &cam);
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

void vcn_mesh_save_png(const vcn_mesh_t *const mesh,
		       const char* filename,
		       int width, int height)
{
	nb_drawing_export_png(filename, width, height, draw_mesh, mesh);
}

void vcn_dewall_save_png(const vcn_mesh_t *const mesh,
			 const char* filename, int width, int height,
			 uint8_t axe, double alpha, uint32_t N,
			 void *vtx_array)
{
	dewall_data_t data;
	data.mesh = mesh;
	data.axe = axe;
	data.alpha = alpha;
	data.N = N;
	data.vtx_array = vtx_array;
	nb_drawing_export_png(filename, width, height, draw_dewall, &data);
}

static void draw_dewall(void *draw_ptr, int width, int height,
			const void *const data_ptr)
{
	const dewall_data_t *const data = data_ptr;

	msh_vtx_t **vertices = data->vtx_array;

	double box[4];
	vcn_utils2D_get_enveloping_box(data->mesh->N_input_vtx,
				       data->mesh->input_vtx,
				       sizeof(*(data->mesh->input_vtx)),
				       msh_vtx_get_x,
				       msh_vtx_get_y,
				       box);

	camera_t cam;
	nb_drawing_utils_set_center_and_zoom(&cam, box, width, height);

	draw_mesh(draw_ptr, width, height, data->mesh);

	/* Draw Wall */
	nb_drawing_set_source_rgb(draw_ptr, 1.0, 0.0, 0.0);
	if (0 == data->axe) {
		double x_alpha = cam.zoom * data->alpha -
			cam.zoom * cam.center[0] + width / 2.0;
		nb_drawing_raw_move_to(draw_ptr, x_alpha, 0);
		nb_drawing_raw_line_to(draw_ptr, x_alpha, height);
	} else {
		double y_alpha = -cam.zoom * data->alpha +
			cam.zoom * cam.center[1] + height / 2.0;
		nb_drawing_raw_move_to(draw_ptr, 0, y_alpha);
		nb_drawing_raw_line_to(draw_ptr, width, y_alpha);
	}
	nb_drawing_stroke(draw_ptr);

	/* Draw vertices of partition */
	nb_drawing_set_source_rgb(draw_ptr, 0.0, 1.0, 0.0);
	nb_drawing_set_line_width(draw_ptr, 1.0);
	draw_vertices_halos(draw_ptr, data->N, vertices, width, height, &cam);
}

static void draw_vertices_halos(void *draw_ptr,
				uint32_t N, msh_vtx_t **vertices,
				int width, int height,
				const camera_t *const cam)
{
	for (uint32_t i = 0; i < N; i++) {
		nb_drawing_move_to(draw_ptr, cam,
				   vertices[i]->x[0],
				   vertices[i]->x[1]);
		nb_drawing_set_circle(draw_ptr, cam,
				      vertices[i]->x[0],
				      vertices[i]->x[1],
				      5.0, true);
		nb_drawing_stroke(draw_ptr);
	}
}

void vcn_mesh_save_eps(const vcn_mesh_t *const mesh,
		   const char* filename,
		   int width, int height)
{
	nb_drawing_export_eps(filename, width, height, draw_mesh, mesh);
}
