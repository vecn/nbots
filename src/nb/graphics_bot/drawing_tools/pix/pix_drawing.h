#ifndef __NB_GRAPHICS_BOT_DRAWING_TOOLS_PIX_PIX_DRAWING_H__
#define __NB_GRAPHICS_BOT_DRAWING_TOOLS_PIX_PIX_DRAWING_H__

#include <stdbool.h>
#include <stdint.h>

#include "nb/graphics_bot/drawing_tools.h"

void* nb_graphics_pix_create_context(int width, int height);
void nb_graphics_pix_destroy_context(void *ctx);
void nb_graphics_pix_export_context(void *ctx, const char *filename);

void nb_graphics_pix_move_to(void *ctx, double x, double y);

void nb_graphics_pix_line_to(void *ctx, double x, double y);

void nb_graphics_pix_arc(void *ctx,
			 double x, double y, double r,
			 double a0, double a1);

void nb_graphics_pix_set_circle(void *ctx,
				double x, double y, double r);

void nb_graphics_pix_set_rectangle(void *ctx, double x1, double y1,
				   double x2, double y2);

void nb_graphics_pix_close_path(void *ctx);

void nb_graphics_pix_set_line_width(void *ctx, double w);

void nb_graphics_pix_set_source_rgb(void *ctx,
				    uint8_t r, uint8_t g, uint8_t b);

void nb_graphics_pix_set_source_rgba(void *ctx, uint8_t r, uint8_t g,
				     uint8_t b, uint8_t a);

void nb_graphics_pix_set_source_grad(void *ctx,
				     nb_graphics_grad_t grad,
				     double x1, double y1,
				     double x2, double y2,
				     nb_graphics_palette_t *pat);

void nb_graphics_pix_set_source_trg(void *ctx,
				    double x1, double y1,
				    double x2, double y2,
				    double x3, double y3,
				    uint8_t r1, uint8_t g1, uint8_t b1,
				    uint8_t r2, uint8_t g2, uint8_t b2,
				    uint8_t r3, uint8_t g3, uint8_t b3);

void nb_graphics_pix_fill(void *ctx);

void nb_graphics_pix_fill_preserve(void *ctx);

void nb_graphics_pix_stroke(void *ctx);

void nb_graphics_pix_set_font_type(void *ctx, const char *type);

void nb_graphics_pix_set_font_size(void *ctx, uint16_t size);

void nb_graphics_pix_show_text(void *ctx, const char *str);
void nb_graphics_pix_get_text_attr(void *ctx, const char *label,
				   nb_text_attr_t *attr);

#endif
