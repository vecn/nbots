#ifndef __NB_GRAPHICS_BOT_DRAWING_TOOLS_H__
#define __NB_GRAPHICS_BOT_DRAWING_TOOLS_H__

#include <stdbool.h>
#include <stdint.h>

#include "nb/graphics_bot/drawing_utils.h"

typedef nb_graphics_context_s nb_graphics_context_t;
typedef void nb_pattern_t;

typedef enum {
	LINEAR, RADIAL
} nb_graphics_grad_t;

typedef struct {
	double width;
	double height;
	double x_left;
	double y_top;
} nb_text_attr_t;

typedef struct {
	int width;
	int height;
	double center[2];
	double zoom;
} camera_t;

void nb_graphics_export_png(const char* filename, int width, int height,
			   void (*draw)(nb_graphics_context_t *g, int w, int h,
					const void *const data),
			   const void *const data);

void nb_graphics_export_eps(const char* filename, int width, int height,
			   void (*draw)(nb_graphics_context_t *g, int w, int h,
					const void *const data),
			   const void *const data);

void nb_graphics_export_asy(const char* filename, int width, int height,
			   void (*draw)(nb_graphics_context_t *g, int w, int h,
					const void *const data),
			   const void *const data);

void nb_graphics_set_camera(nb_graphics_context_t *g, const camera_t *cam);
void nb_graphics_unset_camera(nb_graphics_context_t *g);

void nb_graphics_move_to(nb_graphics_context_t *g, double x, double y);
void nb_graphics_line_to(nb_graphics_context_t *g, double x, double y);
void nb_graphics_arc(nb_graphics_context_t *g, double x, double y, double r,
		    double a0, double a1);

void nb_graphics_set_circle(nb_graphics_context_t *g,
			    double x, double y, double r);

void nb_graphics_set_point(nb_graphics_context_t *g,
			   double x, double y, double size);

void nb_graphics_set_rectangle(nb_graphics_context_t *g, double x1, double y1,
				  double x2, double y2);

void nb_graphics_close_path(nb_graphics_context_t *g);

void nb_graphics_set_line_width(nb_graphics_context_t *g, double w);

void nb_graphics_set_source_rgb(nb_graphics_context_t *g,
				uint8_t r, uint8_t g, uint8_t b);

void nb_graphics_set_source_rgba(nb_graphics_context_t *g,
				 uint8_t r, uint8_t g, uint8_t b, uint8_t a);

void nb_graphics_set_source_grad(nb_graphics_context_t *g,
				 nb_graphics_grad_t grad,
				 double x1, double y1,
				 double x2, double y2,
				 nb_palette_t *pat);

void nb_graphics_set_source_trg(nb_graphics_context_t *g,
				double x1, double y1,
				double x2, double y2,
				double x3, double y3,
				uint8_t r1, uint8_t g1, uint8_t b1,
				uint8_t r2, uint8_t g2, uint8_t b2,
				uint8_t r3, uint8_t g3, uint8_t b3);

void nb_graphics_fill(nb_graphics_context_t *g);

void nb_graphics_fill_preserve(nb_graphics_context_t *g);

void nb_graphics_stroke(nb_graphics_context_t *g);

void nb_graphics_set_font_type(nb_graphics_context_t *g, const char *type);

void nb_graphics_set_font_size(nb_graphics_context_t *g, uint16_t size);

void nb_graphics_show_text(nb_graphics_context_t *g, const char *str);
void nb_graphics_get_text_attr(nb_graphics_context_t *g, const char *label,
			      nb_text_attr_t *attr);


#endif
