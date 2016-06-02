#ifndef __NB_GRAPHICS_BOT_DRAWING_TOOLS_H__
#define __NB_GRAPHICS_BOT_DRAWING_TOOLS_H__

#include <stdbool.h>
#include <stdint.h>

#include "nb/graphics_bot/drawing_utils.h"

void nb_drawing_export_png(const char* filename, int width, int height,
			   void (*draw)(void *draw_ptr, int w, int h,
					const void *const data),
			   const void *const data);

void nb_drawing_export_eps(const char* filename, int width, int height,
			   void (*draw)(void *draw_ptr, int w, int h,
					const void *const data),
			   const void *const data);

void nb_drawing_raw_move_to(void *draw_ptr, double x, double y);

void nb_drawing_raw_line_to(void *draw_ptr, double x, double y);

void nb_drawing_move_to(void *draw_ptr, const camera_t *cam,
			double x, double y);

void nb_drawing_line_to(void *draw_ptr, const camera_t *cam,
			double x, double y);

void nb_drawing_arc(void *draw_ptr, const camera_t *cam,
		    double x, double y, double r,
		    double a0, double a1);

void nb_drawing_set_circle(void *draw_ptr, const camera_t *cam,
			   double x, double y, double r, bool r_fixed);

void nb_drawing_close_path(void *draw_ptr);

void nb_drawing_set_line_width(void *draw_ptr, double w);

void nb_drawing_set_source_rgb(void *draw_ptr, double r,
			       double g, double b);

void nb_drawing_set_source_rgba(void *draw_ptr, double r,
				double g, double b, double a);

void nb_drawing_fill(void *draw_ptr);

void nb_drawing_fill_preserve(void *draw_ptr);

void nb_drawing_stroke(void *draw_ptr);

void nb_drawing_paint(void *draw_ptr);

void nb_drawing_set_font_type(void *draw_ptr, const char *type);

void nb_drawing_set_font_size(void *draw_ptr, uint16_t size);

void nb_drawing_show_text(void *draw_ptr, const char *str);

#endif
