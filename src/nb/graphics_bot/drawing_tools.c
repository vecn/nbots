/* http://paulbourke.net/dataformats/postscript/ */

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

void nb_drawing_raw_set_rectangle(void *draw_ptr, double x1, double y1,
				  double x2, double y2)
{
	cairo_rectangle(draw_ptr, x1, y1, x2, y2);
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

void nb_drawing_set_source(void *draw_ptr, nb_pattern_t *pat)
{
	cairo_set_source(draw_ptr, pat);
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

void nb_drawing_get_text_attr(void *draw_ptr, const char *label,
			      nb_text_attr_t *attr)
{
	cairo_text_extents_t extents;
	cairo_text_extents(draw_ptr, label, &extents);
	attr->width = extents.width;
	attr->height = extents.height;
	attr->x_left = extents.x_bearing;
	attr->y_top = extents.y_bearing;
}

nb_pattern_t* nb_pattern_create(void)
{
	return cairo_pattern_create_mesh();
}

nb_pattern_t* nb_pattern_create_linear(double x1, double y1,
				       double x2, double y2)
{
	cairo_pattern_create_linear(x1, y1, x2, y2);
}

void nb_pattern_destroy(nb_pattern_t *pat)
{
	cairo_pattern_destroy(pat);
}

void nb_pattern_add_color_stop_rgb(nb_pattern_t *pat, float stop,
				   double r, double g, double b)
{
	cairo_pattern_add_color_stop_rgb(pat, stop, r, g, b);
}

void nb_pattern_begin_patch(nb_pattern_t *pat)
{
	cairo_mesh_pattern_begin_patch(pat);
}

void nb_pattern_move_to(nb_pattern_t *pat, const camera_t *cam,
			double x, double y)
{
	double cam_x = cam->zoom * (x - cam->center[0]) + cam->width / 2.0;
	double cam_y = cam->zoom * (cam->center[1] - y) + cam->height / 2.0;
	cairo_mesh_pattern_move_to(pat, cam_x, cam_y);
}

void nb_pattern_line_to(nb_pattern_t *pat, const camera_t *cam,
			double x, double y)
{
	double cam_x = cam->zoom * (x - cam->center[0]) + cam->width / 2.0;
	double cam_y = cam->zoom * (cam->center[1] - y) + cam->height / 2.0;
	cairo_mesh_pattern_line_to(pat, cam_x, cam_y);
}

void nb_pattern_end_patch(nb_pattern_t *pat)
{
	cairo_mesh_pattern_end_patch(pat);
}

void nb_pattern_set_corner_color_rgb(nb_pattern_t *pat, int vtx_id,
				     double r, double g, double b)
{
	cairo_mesh_pattern_set_corner_color_rgb(pat, vtx_id, r, g, b);
}
