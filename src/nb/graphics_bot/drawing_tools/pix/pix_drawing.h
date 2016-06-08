#ifndef __NB_GRAPHICS_BOT_DRAWING_TOOLS_PIX_PIX_DRAWING_H__
#define __NB_GRAPHICS_BOT_DRAWING_TOOLS_PIX_PIX_DRAWING_H__

#include <stdbool.h>
#include <stdint.h>

#include "nb/graphics_bot/drawing_tools.h"

void nb_graphics_pix_move_to(void *draw_ptr, double x, double y);

void nb_graphics_pix_line_to(void *draw_ptr, double x, double y);

void nb_graphics_pix_arc(void *draw_ptr,
			 double x, double y, double r,
			 double a0, double a1);

void nb_graphics_pix_set_circle(void *draw_ptr,
				double x, double y, double r);

void nb_graphics_pix_set_rectangle(void *draw_ptr, double x1, double y1,
				   double x2, double y2);

void nb_graphics_pix_close_path(void *draw_ptr);

void nb_graphics_pix_set_line_width(void *draw_ptr, double w);

void nb_graphics_pix_set_source_rgb(void *draw_ptr,
				    uint8_t r, uint8_t g, uint8_t b);

void nb_graphics_pix_set_source_rgba(void *draw_ptr, uint8_t r, uint8_t g,
				     uint8_t b, uint8_t a);

void nb_graphics_pix_set_source(void *draw_ptr, nb_palette_t *pat);

void nb_graphics_pix_fill(void *draw_ptr);

void nb_graphics_pix_fill_preserve(void *draw_ptr);

void nb_graphics_pix_stroke(void *draw_ptr);

void nb_graphics_pix_set_font_type(void *draw_ptr, const char *type);

void nb_graphics_pix_set_font_size(void *draw_ptr, uint16_t size);

void nb_graphics_pix_show_text(void *draw_ptr, const char *str);
void nb_graphics_pix_get_text_attr(void *draw_ptr, const char *label,
				   nb_text_attr_t *attr);

#endif
