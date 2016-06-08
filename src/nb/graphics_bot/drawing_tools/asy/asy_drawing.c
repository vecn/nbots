#include <stdint.h>

#include "nb/math_bot.h"
#include "nb/graphics_bot/drawing_tools.h"

#include "asy_drawing.h"

void nb_graphics_asy_move_to(void *draw_ptr, double x, double y)
{
	;
}

void nb_graphics_asy_line_to(void *draw_ptr, double x, double y)
{
	;
}

void nb_graphics_asy_arc(void *draw_ptr,
			 double x, double y, double r,
			 double a0, double a1)
{
	;
}

void nb_graphics_asy_set_circle(void *draw_ptr,
				double x, double y, double r)
{

}

void nb_graphics_asy_set_rectangle(void *draw_ptr, double x1, double y1,
				       double x2, double y2)
{
	;
}

void nb_graphics_asy_close_path(void *draw_ptr)
{
	;
}

void nb_graphics_asy_set_line_width(void *draw_ptr, double w)
{
	;
}

void nb_graphics_asy_set_source_rgb(void *draw_ptr,
				    uint8_t r, uint8_t g, uint8_t b)
{
	;
}

void nb_graphics_asy_set_source_rgba(void *draw_ptr, uint8_t r, uint8_t g,
				     uint8_t b, uint8_t a)
{
	;
}

void nb_graphics_asy_set_source(void *draw_ptr, nb_palette_t *pat)
{
	;
}

void nb_graphics_asy_fill(void *draw_ptr)
{
	;
}

void nb_graphics_asy_fill_preserve(void *draw_ptr)
{
	;
}

void nb_graphics_asy_stroke(void *draw_ptr)
{
	;
}

void nb_graphics_asy_set_font_type(void *draw_ptr, const char *type)
{
	;
}

void nb_graphics_asy_set_font_size(void *draw_ptr, uint16_t size)
{
	;
}

void nb_graphics_asy_show_text(void *draw_ptr, const char *str)
{
	;
}

void nb_graphics_asy_get_text_attr(void *draw_ptr, const char *label,
				   nb_text_attr_t *attr)
{
	;
}
