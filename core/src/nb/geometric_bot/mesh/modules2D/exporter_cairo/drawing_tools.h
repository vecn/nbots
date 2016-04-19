#ifndef __NB_GEOMETRIC_BOT_MESH_MODULES2D_EXPORTER_CAIRO_DRAWING_TOOLS_H__
#define __NB_GEOMETRIC_BOT_MESH_MODULES2D_EXPORTER_CAIRO_DRAWING_TOOLS_H__

#include "drawing_utils.h"

void nb_drawing_export_png(const char* filename, int width, int height,
			   void (*draw)(void *draw_ptr, int w, int h,
					const void *const data),
			   const void *const data);

void nb_drawing_export_eps(const char* filename, int width, int height,
			   void (*draw)(void *draw_ptr, int w, int h,
					const void *const data),
			   const void *const data);

void nb_drawing_move_to(void *draw_ptr, const camera_t *cam,
			double x, double y);

void nb_drawing_line_to(void *draw_ptr, const camera_t *cam,
			double x, double y);

void nb_drawing_arc(void *draw_ptr, const camera_t *cam,
		    double x, double y, double r,
		    double a0, double a1);

void nb_drawing_close_path(void *draw_ptr);

void nb_drawing_set_line_width(void *draw_ptr, double w);

void nb_drawing_set_source_rgb(void *draw_ptr, double r,
			       double g, double b);

void nb_drawing_set_source_rgba(void *draw_ptr, double r,
				double g, double b, double a);

void nb_drawing_fill(void *draw_ptr);
void nb_drawing_fill_preserve(void *draw_ptr);
void nb_drawing_stroke(void *draw_ptr);

#endif
