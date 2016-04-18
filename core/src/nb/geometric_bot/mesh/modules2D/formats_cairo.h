
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
				 const void *const data);
