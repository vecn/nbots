#include <stdint.h>

#include "nb/math_bot.h"
#include "nb/graphics_bot/drawing_tools.h"

#include "eps_drawing.h"


void* nb_graphics_eps_create_context(int width, int height)
{
	return 0;
}

void nb_graphics_eps_destroy_context(void *ctx)
{
	;
}

void nb_graphics_eps_export_context(const void *ctx, const char *filename)
{
	;
}

void nb_graphics_eps_move_to(void *ctx, float x, float y)
{
	;
}

void nb_graphics_eps_line_to(void *ctx, float x, float y)
{
	;
}

void nb_graphics_eps_arc(void *ctx,
			 float x, float y, float r,
			 float a0, float a1)
{
	;
}

void nb_graphics_eps_close_path(void *ctx)
{
	;
}

void nb_graphics_eps_set_circle(void *ctx,
				float x, float y, float r)
{

}

void nb_graphics_eps_set_rectangle(void *ctx, float x1, float y1,
				       float x2, float y2)
{
	;
}

void nb_graphics_eps_set_line_width(void *ctx, float w)
{
	;
}

void nb_graphics_eps_set_source_rgb(void *ctx,
				    uint8_t r, uint8_t g, uint8_t b)
{
	;
}

void nb_graphics_eps_set_source_rgba(void *ctx, uint8_t r, uint8_t g,
				     uint8_t b, uint8_t a)
{
	;
}

void nb_graphics_eps_set_source_grad(void *ctx,
				     nb_graphics_grad_t grad,
				     float x1, float y1,
				     float x2, float y2,
				     nb_graphics_palette_t *pat)
{
	;
}

void nb_graphics_eps_set_source_trg(void *ctx,
				    float x1, float y1,
				    float x2, float y2,
				    float x3, float y3,
				    const uint8_t rgba1[4],
				    const uint8_t rgba2[4],
				    const uint8_t rgba3[4])
{
	;
}

void nb_graphics_eps_fill(void *ctx)
{
	;
}

void nb_graphics_eps_fill_preserve(void *ctx)
{
	;
}

void nb_graphics_eps_stroke(void *ctx)
{
	;
}

void nb_graphics_eps_set_font_type(void *ctx, const char *type)
{
	;
}

void nb_graphics_eps_set_font_size(void *ctx, uint16_t size)
{
	;
}

void nb_graphics_eps_show_text(void *ctx, const char *str)
{
	;
}

void nb_graphics_eps_get_text_attr(const void *ctx, const char *label,
				   nb_text_attr_t *attr)
{
	;
}
