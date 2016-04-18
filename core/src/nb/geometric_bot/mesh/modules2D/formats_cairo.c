#include <cairo.h>
#include <cairo-ps.h>

#include "formats_cairo.h"

void nb_cairo_drawing_export_png(const char* filename, int width, int height,
				 void (*draw)(cairo_t *cr, int w, int h,
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

void nb_cairo_drawing_export_eps(const char* filename, int width, int height,
				 void (*draw)(cairo_t *cr, int w, int h,
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
