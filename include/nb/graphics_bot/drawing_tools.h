#ifndef __NB_GRAPHICS_BOT_DRAWING_TOOLS_H__
#define __NB_GRAPHICS_BOT_DRAWING_TOOLS_H__

#include <stdbool.h>
#include <stdint.h>

#include "nb/graphics_bot/drawing_utils.h"

typedef void nb_pattern_t;
typedef struct {
	double width;
	double height;
	double x_left;
	double y_top;
} nb_text_attr_t;

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

void nb_drawing_raw_set_rectangle(void *draw_ptr, double x1, double y1,
				  double x2, double y2);

void nb_drawing_close_path(void *draw_ptr);

void nb_drawing_set_line_width(void *draw_ptr, double w);

void nb_drawing_set_source_rgb(void *draw_ptr, double r,
			       double g, double b);

void nb_drawing_set_source_rgba(void *draw_ptr, double r,
				double g, double b, double a);

void nb_drawing_set_source(void *draw_ptr, nb_pattern_t *pat);

void nb_drawing_fill(void *draw_ptr);

void nb_drawing_fill_preserve(void *draw_ptr);

void nb_drawing_stroke(void *draw_ptr);

void nb_drawing_paint(void *draw_ptr);

void nb_drawing_set_font_type(void *draw_ptr, const char *type);

void nb_drawing_set_font_size(void *draw_ptr, uint16_t size);

void nb_drawing_show_text(void *draw_ptr, const char *str);
void nb_drawing_get_text_attr(void *draw_ptr, const char *label,
			      nb_text_attr_t *attr);

nb_pattern_t* nb_pattern_create(void);
nb_pattern_t* nb_pattern_create_linear(double x1, double y1,
				       double x2, double y2);
void nb_pattern_destroy(nb_pattern_t *pat);
void nb_pattern_add_color_stop_rgb(nb_pattern_t *pat, float stop,
				   double r, double g, double b);
void nb_pattern_begin_patch(nb_pattern_t *pat);
void nb_pattern_move_to(nb_pattern_t *pat, const camera_t *cam,
			double x, double y);
void nb_pattern_line_to(nb_pattern_t *pat, const camera_t *cam,
			double x, double y);
void nb_pattern_end_patch(nb_pattern_t *pat);
void nb_pattern_set_corner_color_rgb(nb_pattern_t *pat, int vtx_id,
				     double r, double g, double b);

#endif
