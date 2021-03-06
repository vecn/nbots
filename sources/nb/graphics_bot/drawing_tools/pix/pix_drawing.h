#ifndef __NB_GRAPHICS_BOT_DRAWING_TOOLS_PIX_PIX_DRAWING_H__
#define __NB_GRAPHICS_BOT_DRAWING_TOOLS_PIX_PIX_DRAWING_H__

#include <stdbool.h>
#include <stdint.h>

#include "nb/graphics_bot/drawing_tools.h"

void* nb_graphics_pix_create_context(int width, int height);
void nb_graphics_pix_destroy_context(void *ctx);
void nb_graphics_pix_export_context(const void *ctx, const char *filename);

void nb_graphics_pix_move_to(void *ctx, float x, float y);

void nb_graphics_pix_line_to(void *ctx, float x, float y);

void nb_graphics_pix_qcurve_to(void *ctx, float x, float y,
			       float xcontrol, float ycontrol);
void nb_graphics_pix_qrcurve_to(void *ctx, float x, float y,
				float xcontrol, float ycontrol, float w);
void nb_graphics_pix_curve_to(void *ctx, float x, float y,
			      float x0_control, float y0_control,
			      float x1_control, float y1_control);
void nb_graphics_pix_close_path(void *ctx);

void nb_graphics_pix_set_line_width(void *ctx, float w);

void nb_graphics_pix_set_source_rgb(void *ctx,
				    uint8_t r, uint8_t g, uint8_t b);

void nb_graphics_pix_set_source_rgba(void *ctx, uint8_t r, uint8_t g,
				     uint8_t b, uint8_t a);

void nb_graphics_pix_set_source_grad(void *ctx,
				     nb_graphics_grad_t grad,
				     float x1, float y1,
				     float x2, float y2,
				     nb_palette_t *pat);

void nb_graphics_pix_set_source_trg(void *ctx,
				    float x1, float y1,
				    float x2, float y2,
				    float x3, float y3,
				    const uint8_t rgba1[4],
				    const uint8_t rgba2[4],
				    const uint8_t rgba3[4]);

void nb_graphics_pix_fill(void *ctx);

void nb_graphics_pix_fill_preserve(void *ctx);

void nb_graphics_pix_stroke(void *ctx);

void nb_graphics_pix_stroke_preserve(void *ctx);

void nb_graphics_pix_set_font_type(void *ctx, const char *type);

void nb_graphics_pix_set_font_size(void *ctx, uint16_t size);

void nb_graphics_pix_show_text(void *ctx, int x, int y, const char *str);
void nb_graphics_pix_get_text_attr(const void *ctx, const char *label,
				   nb_graphics_text_attr_t *attr);

#endif
