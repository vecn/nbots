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

static void draw_triangles(void *const g,
			   const nb_container_t *const ht_trg);
static void draw_input_segments(void *const g,
				uint32_t N_input_segments,
				msh_edge_t **segments);
static void draw_mesh(nb_graphics_context_t *g, int width, int height,
		      const void *const mesh_ptr);
static double msh_vtx_get_x(const void *const vtx_ptr);
static double msh_vtx_get_y(const void *const vtx_ptr);

static void draw_dewall(nb_graphics_context_t *g, int width, int height,
			const void *const mesh_ptr);
static void draw_vertices_halos(nb_graphics_context_t *g, 
				uint32_t N, msh_vtx_t **vertices);


static void draw_triangles(void *const g,
			   const nb_container_t *const ht_trg)
{
	nb_graphics_set_line_width(g, 0.5);
	nb_iterator_t* iter = nb_iterator_create();
	nb_iterator_set_container(iter, ht_trg);
	while (nb_iterator_has_more(iter)) {
		const msh_trg_t* trg = nb_iterator_get_next(iter);
		nb_graphics_move_to(g, trg->v1->x[0], trg->v1->x[1]);
		nb_graphics_line_to(g, trg->v2->x[0], trg->v2->x[1]);
		nb_graphics_line_to(g, trg->v3->x[0], trg->v3->x[1]);
		nb_graphics_close_path(g);
		nb_graphics_set_source_rgba(g, 25, 75, 255, 128);
		nb_graphics_fill_preserve(g);
		nb_graphics_set_source(g, NB_BLUE);
		nb_graphics_stroke(g);
	}
	nb_iterator_destroy(iter);
}

static void draw_input_segments(void *const g,
				uint32_t N_input_segments,
				msh_edge_t **segments)
{
	nb_graphics_set_source(g, NB_AQUAMARIN);
	nb_graphics_set_line_width(g, 0.75);
	for (uint32_t i = 0; i < N_input_segments; i++) {
		if (NULL == segments[i])
			continue;

		msh_edge_t* sgm = segments[i];

		nb_graphics_move_to(g, sgm->v1->x[0], sgm->v1->x[1]);
		msh_vtx_t* vtx = sgm->v2;
		while (NULL != sgm) {
			nb_graphics_line_to(g, vtx->x[0], vtx->x[1]);
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
		nb_graphics_stroke(g);
	}
}

static void draw_mesh(nb_graphics_context_t *g, int width, int height,
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

	nb_graphics_enable_camera(g);
	camera_t* cam = nb_graphics_get_camera(g);
	nb_graphics_cam_fit_box(cam, box, width, height);

	draw_triangles(g, mesh->ht_trg);

	draw_input_segments(g, mesh->N_input_sgm, mesh->input_sgm);
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

void vcn_mesh_draw(const vcn_mesh_t *const mesh,
		   const char* filename,
		   int width, int height)
{
	nb_graphics_export(filename, width, height, draw_mesh, mesh);
}

void vcn_dewall_draw(const vcn_mesh_t *const mesh,
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
	nb_graphics_export(filename, width, height, draw_dewall, &data);
}

static void draw_dewall(nb_graphics_context_t *g, int width, int height,
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

	nb_graphics_enable_camera(g);
	camera_t* cam = nb_graphics_get_camera(g);
	nb_graphics_cam_fit_box(cam, box, width, height);

	draw_mesh(g, width, height, data->mesh);

	/* Draw Wall */
	nb_graphics_set_source(g, NB_RED);
	if (0 == data->axe) {
		double x_alpha = cam->zoom * data->alpha -
			cam->zoom * cam->center[0] + width / 2.0;
		nb_graphics_move_to(g, x_alpha, 0);
		nb_graphics_line_to(g, x_alpha, height);
	} else {
		double y_alpha = -cam->zoom * data->alpha +
			cam->zoom * cam->center[1] + height / 2.0;
		nb_graphics_move_to(g, 0, y_alpha);
		nb_graphics_line_to(g, width, y_alpha);
	}
	nb_graphics_stroke(g);

	draw_vertices_halos(g, data->N, vertices);
}

static void draw_vertices_halos(nb_graphics_context_t *g,
				uint32_t N, msh_vtx_t **vertices)
{
	nb_graphics_set_source(g, NB_GREEN);
	nb_graphics_set_line_width(g, 1.0);
	for (uint32_t i = 0; i < N; i++) {
		nb_graphics_move_to(g,
				    vertices[i]->x[0],
				    vertices[i]->x[1]);
		nb_graphics_set_point(g,
				      vertices[i]->x[0],
				      vertices[i]->x[1],
				      10.0);
		nb_graphics_stroke(g);
	}
}
