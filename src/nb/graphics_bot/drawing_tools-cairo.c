#include <stdint.h>
#include <cairo.h>
#include <cairo-ps.h>

#include "nb/math_bot.h"
#include "nb/graphics_bot/drawing_tools.h"

void nb_drawing_export_png(const char* filename, int width, int height,
			   void (*draw)(void *draw_ptr, int w, int h,
					const void *const data),
			   const void *const data)
{
	cairo_surface_t* surface =
		cairo_image_surface_create(CAIRO_FORMAT_ARGB32, 
					   width, height);
	cairo_t* cr = cairo_create(surface);

	draw(cr, width, height, data);

	cairo_surface_write_to_png(surface, filename);

	cairo_destroy(cr);
	cairo_surface_destroy(surface);
}

void nb_drawing_export_eps(const char* filename, int width, int height,
			   void (*draw)(void *draw_ptr, int w, int h,
					const void *const data),
			   const void *const data)
{
	/* Create drawable surface and cairo context */
	cairo_surface_t* surface = 
		cairo_ps_surface_create (filename, width, height);
	cairo_ps_surface_set_eps(surface, 1 /* TRUE from cairo_bool_t*/);

	cairo_t* cr = cairo_create(surface);

	/* Initialize Post script */
	cairo_ps_surface_dsc_begin_page_setup(surface);

	if(height < width)
		cairo_ps_surface_dsc_comment(surface,
					     "%%PageOrientation: Portrait");
	else
		cairo_ps_surface_dsc_comment(surface,
					     "%%PageOrientation: Landscape");

	draw(cr, width, height, data);

	cairo_surface_show_page(surface);/* Show post-script */

	cairo_destroy(cr);
	cairo_surface_finish(surface);
	cairo_surface_destroy(surface);
}

void nb_drawing_raw_move_to(void *draw_ptr, double x, double y)
{
	cairo_move_to(draw_ptr, x, y);
}

void nb_drawing_raw_line_to(void *draw_ptr, double x, double y)
{
	cairo_line_to(draw_ptr, x, y);
}

void nb_drawing_move_to(void *draw_ptr, const camera_t *cam,
			double x, double y)
{
	double cam_x = cam->zoom * (x - cam->center[0]) + cam->width / 2.0;
	double cam_y = cam->zoom * (cam->center[1] - y) + cam->height / 2.0;
	cairo_move_to(draw_ptr, cam_x, cam_y);
}

void nb_drawing_line_to(void *draw_ptr, const camera_t *cam,
			double x, double y)
{
	double cam_x = cam->zoom * (x - cam->center[0]) + cam->width / 2.0;
	double cam_y = cam->zoom * (cam->center[1] - y) + cam->height / 2.0;
	cairo_line_to(draw_ptr, cam_x, cam_y);
}

void nb_drawing_arc(void *draw_ptr, const camera_t *cam,
		    double x, double y, double r,
		    double a0, double a1)
{
	double cam_x = cam->zoom * (x - cam->center[0]) + cam->width / 2.0;
	double cam_y = cam->zoom * (cam->center[1] - y) + cam->height / 2.0;
	cairo_arc(draw_ptr, cam_x, cam_y, r, a0, a1);
}

void nb_drawing_set_circle(void *draw_ptr, const camera_t *cam,
			   double x, double y, double r, bool r_fixed)
{
	double cam_x = cam->zoom * (x - cam->center[0]) + cam->width / 2.0;
	double cam_y = cam->zoom * (cam->center[1] - y) + cam->height / 2.0;
	if (!r_fixed)
		r *= cam->zoom;

	cairo_move_to(draw_ptr, cam_x + r, cam_y);
	cairo_arc(draw_ptr, cam_x, cam_y, r, 0, 2.0 * NB_PI);
}

void nb_drawing_close_path(void *draw_ptr)
{
	cairo_close_path(draw_ptr);
}

void nb_drawing_set_line_width(void *draw_ptr, double w)
{
	cairo_set_line_width(draw_ptr, w);
}

void nb_drawing_set_source_rgb(void *draw_ptr, double r,
			       double g, double b)
{
	cairo_set_source_rgb(draw_ptr, r, g, b);
}

void nb_drawing_set_source_rgba(void *draw_ptr, double r,
				double g, double b, double a)
{
	cairo_set_source_rgba(draw_ptr, r, g, b, a);
}

void nb_drawing_fill(void *draw_ptr)
{
	cairo_fill(draw_ptr);
}

void nb_drawing_fill_preserve(void *draw_ptr)
{
	cairo_fill_preserve(draw_ptr);
}

void nb_drawing_stroke(void *draw_ptr)
{
	cairo_stroke(draw_ptr);
}

void nb_drawing_paint(void *draw_ptr)
{
	cairo_paint(draw_ptr);
}

void nb_drawing_set_font_type(void *draw_ptr, const char *type)
{
	cairo_select_font_face(draw_ptr, type,
			       CAIRO_FONT_SLANT_NORMAL,
			       CAIRO_FONT_WEIGHT_NORMAL);
}

void nb_drawing_set_font_size(void *draw_ptr, uint16_t size)
{
	cairo_set_font_size(draw_ptr, size);
}

void nb_drawing_show_text(void *draw_ptr, const char *str)
{
	cairo_show_text(draw_ptr, str);
}
