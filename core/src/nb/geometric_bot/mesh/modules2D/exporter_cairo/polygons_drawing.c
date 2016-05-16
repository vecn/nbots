#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>

#include "nb/math_bot.h"
#include "nb/container_bot.h"
#include "nb/geometric_bot.h"

#include "drawing_utils.h"
#include "drawing_tools.h"

#include "../../mesh2D_structs.h"

static void draw_mesh(void *draw_ptr, int width, int height,
		      const void *const poly_ptr);
static void draw_polygons(void *draw_ptr, const nb_mshpoly_t *const poly,
			  const camera_t *cam);
static void draw_centroids(void *draw_ptr, const nb_mshpoly_t *const poly,
			  const camera_t *cam);
static void draw_input_sgm(void *draw_ptr, const nb_mshpoly_t *const poly,
			  const camera_t *cam);

void nb_mshpoly_export_png(const nb_mshpoly_t *const poly,
			   const char* filename, int width, int height)
{
	nb_drawing_export_png(filename, width, height, draw_mesh, poly);
}

void nb_mshpoly_export_eps(const nb_mshpoly_t *const poly,
			   const char* filename, int width, int height)
{
	nb_drawing_export_eps(filename, width, height, draw_mesh, poly);
}


static void draw_mesh(void *draw_ptr, int width, int height,
		      const void *const poly_ptr)
{
	const nb_mshpoly_t *const poly = poly_ptr;

	double box[4];
	vcn_utils2D_get_enveloping_box(poly->N_nod, poly->nod,
				       2 * sizeof(*(poly->nod)),
				       vcn_utils2D_get_x_from_darray,
				       vcn_utils2D_get_y_from_darray,
				       box);
	camera_t cam;
	nb_drawing_utils_set_center_and_zoom(&cam, box, width, height);
	draw_polygons(draw_ptr, poly, &cam);
	draw_centroids(draw_ptr, poly, &cam);
	draw_input_sgm(draw_ptr, poly, &cam);
}

static void draw_polygons(void *draw_ptr, const nb_mshpoly_t *const poly,
			  const camera_t *cam)
{
	nb_drawing_set_line_width(draw_ptr, 0.5);
	for (uint32_t i = 0; i < poly->N_elems; i++) {
		uint32_t id = poly->adj[i][0];
		double x = poly->nod[id * 2];
		double y = poly->nod[id*2+1];
		nb_drawing_move_to(draw_ptr, cam, x, y);
		for (uint16_t j = 1; j < poly->N_adj; j++) {
			id = poly->adj[i][j];
			x = poly->nod[id * 2];
			y = poly->nod[id*2+1];
			nb_drawing_line_to(draw_ptr, cam, x, y);			
		}
		nb_drawing_close_path(draw_ptr);

		nb_drawing_set_source_rgba(draw_ptr, 0.2, 0.9, 0.9, 0.5);
		nb_drawing_fill_preserve(draw_ptr);

		nb_drawing_set_source_rgb(draw_ptr, 0, 0, 0);
		nb_drawing_stroke(draw_ptr);
	}
}

static void draw_centroids(void *draw_ptr, const nb_mshpoly_t *const poly,
			  const camera_t *cam)
{
	double r = 3.0;
	nb_drawing_set_source_rgb(draw_ptr, 0.2, 0.2, 9.0);
	for (uint32_t i = 0; i < poly->N_elems; i++) {;
		double x = poly->cen[i * 2];
		double y = poly->cen[i*2+1];
		nb_drawing_move_to(draw_ptr, cam, x + r, y);
		nb_drawing_arc(draw_ptr, cam, x, y, r, 0, 2.0 * NB_PI);

		nb_drawing_fill(draw_ptr);
	}
	
}

static void draw_input_sgm(void *draw_ptr, const nb_mshpoly_t *const poly,
			  const camera_t *cam)
{
	nb_drawing_set_line_width(draw_ptr, 1.0);
	nb_drawing_set_source_rgb(draw_ptr, 0.5, 0.3, 1.0);
	for (uint32_t i = 0; i < poly->N_sgm; i++) {
		if (0 < poly->N_nod_x_sgm[i]) {
			uint32_t nj = poly->nod_x_sgm[i][0];
			nb_drawing_move_to(draw_ptr, cam,
					   poly->nod[nj * 2],
					   poly->nod[nj*2+1]);
			for (uint32_t j = 1; j < poly->N_nod_x_sgm[i]; j++) {
				nj = poly->nod_x_sgm[i][j];
				nb_drawing_line_to(draw_ptr, cam,
						   poly->nod[nj * 2],
						   poly->nod[nj*2+1]);
			}
			nb_drawing_stroke(draw_ptr);
		}
	}

}
