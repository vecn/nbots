#include <stdio.h>
#include <stdint.h>
#include <alloca.h>

#include "nb/math_bot.h"
#include "nb/graphics_bot/drawing_tools.h"

#include "drawing_tools/pix/pix_drawing.h";
#include "drawing_tools/eps/eps_drawing.h";
#include "drawing_tools/asy/asy_drawing.h";

static void set_pix_tools(nb_graphics_context_t *g);
static void set_eps_tools(nb_graphics_context_t *g);
static void set_asy_tools(nb_graphics_context_t *g);
static void get_camera_view(const camera_t *cam, double cam_vtx[2],
			    double x, double y);

struct nb_graphics_context_s {
	const camera_t *cam;
	void *draw_ptr;
	void (*move_to)(void *draw_ptr, double x, double y);
	void (*line_to)(void *draw_ptr, double x, double y);
	void (*move_to)(void *draw_ptr, double x, double y);
	void (*line_to)(void *draw_ptr, double x, double y);
	void (*arc)(void *draw_ptr, double x, double y, double r,
		    double a0, double a1);
	void (*set_circle)(void *draw_ptr, double x, double y, double r);
	void (*set_rectangle)(void *draw_ptr, double x1, double y1,
				  double x2, double y2);
	void (*close_path)(void *draw_ptr);
	void (*set_line_width)(void *draw_ptr, double w);
	void (*set_source_rgb)(void *draw_ptr, uint8_t r, uint8_t g, uint8_t b);
	void (*set_source_rgba)(void *draw_ptr, uint8_t r, uint8_t g,
				uint8_t b, uint8_t a);
	void (*set_source_grad)(nb_graphics_context_t *g,
				nb_graphics_grad_t grad,
				double x1, double y1,
				double x2, double y2,
				nb_palette_t *pat);
	void (*set_source_trg)(nb_graphics_context_t *g,
			       double x1, double y1,
			       double x2, double y2,
			       double x3, double y3,
			       uint8_t r1, uint8_t g1, uint8_t b1,
			       uint8_t r2, uint8_t g2, uint8_t b2,
			       uint8_t r3, uint8_t g3, uint8_t b3);
	void (*fill)(void *draw_ptr);
	void (*fill_preserve)(void *draw_ptr);
	void (*stroke)(void *draw_ptr);
	void (*set_font_type)(void *draw_ptr, const char *type);
	void (*set_font_size)(void *draw_ptr, uint16_t size);
	void (*show_text)(void *draw_ptr, const char *str);
	void (*get_text_attr)(void *draw_ptr, const char *label,
			      nb_text_attr_t *attr);
};

void nb_graphics_export_png(const char* filename, int width, int height,
			   void (*draw)(nb_graphics_context_t *g, int w, int h,
					const void *const data),
			   const void *const data)
{
	uint32_t memsize = vcn_image_get_memsize();
	vcn_image_t *img = alloca(memsize);
	vcn_image_init_white(img, width, height, 4);
	
	nb_graphics_context_t g;
	set_pix_tools(&g);
	g.cam = NULL;
	g.draw_ptr = img;
	draw(&g, width, height, data);

	vcn_image_write_png(img);
	vcn_image_finish(img);
	fclose(fp);
}

static void set_pix_tools(nb_graphics_context_t *g)
{
	g->move_to = nb_graphics_pix_move_to;
	g->line_to = nb_graphics_pix_line_to;
	g->arc = nb_graphics_pix_arc;
	g->set_circle = nb_graphics_pix_set_circle;
	g->set_rectangle = nb_graphics_pix_set_rectangle;
	g->close_path = nb_graphics_pix_close_path;
	g->set_line_width = nb_graphics_pix_line_width;
	g->set_source_rgb = nb_graphics_pix_source_rgb;
	g->set_source_rgba = nb_graphics_pix_source_rgba;
	g->set_source_grad = nb_graphics_pix_source_grad;
	g->set_source_trg = nb_graphics_pix_source_trg;
	g->fill = nb_graphics_pix_fill;
	g->fill_preserve = nb_graphics_pix_fill_preserve;
	g->stroke = nb_graphics_pix_stroke;
	g->set_font_type = nb_graphics_pix_set_font_type;
	g->set_font_size = nb_graphics_pix_set_font_size;
	g->show_text = nb_graphics_pix_show_text;
	g->get_text_attr = nb_graphics_pix_get_text_attr;
}

void nb_graphics_export_eps(const char* filename, int width, int height,
			   void (*draw)(nb_graphics_context_t *g, int w, int h,
					const void *const data),
			   const void *const data)
{
	FILE *fp = fopen(filename, "w");
	if (NULL != fp) {
		nb_graphics_context_t g;
		set_eps_tools(&g);
		g.cam = NULL;
		g.draw_ptr = fp;
		draw(&g, width, height, data);
	}
	fclose(fp);
}

static void set_eps_tools(nb_graphics_context_t *g)
{
	g->move_to = nb_graphics_eps_move_to;
	g->line_to = nb_graphics_eps_line_to;
	g->arc = nb_graphics_eps_arc;
	g->set_circle = nb_graphics_eps_set_circle;
	g->set_rectangle = nb_graphics_eps_set_rectangle;
	g->close_path = nb_graphics_eps_close_path;
	g->set_line_width = nb_graphics_eps_line_width;
	g->set_source_rgb = nb_graphics_eps_source_rgb;
	g->set_source_rgba = nb_graphics_eps_source_rgba;
	g->set_source_grad = nb_graphics_eps_source_grad;
	g->set_source_trg = nb_graphics_eps_source_trg;
	g->fill = nb_graphics_eps_fill;
	g->fill_preserve = nb_graphics_eps_fill_preserve;
	g->stroke = nb_graphics_eps_stroke;
	g->set_font_type = nb_graphics_eps_set_font_type;
	g->set_font_size = nb_graphics_eps_set_font_size;
	g->show_text = nb_graphics_eps_show_text;
	g->get_text_attr = nb_graphics_eps_get_text_attr;
}

void nb_graphics_export_asy(const char* filename, int width, int height,
			   void (*draw)(nb_graphics_context_t *g, int w, int h,
					const void *const data),
			   const void *const data)
{
	FILE *fp = fopen(filename, "w");
	if (NULL != fp) {
		nb_graphics_context_t g;
		set_asy_tools(&g);
		g.cam = NULL;
		g.draw_ptr = fp;
		draw(&g, width, height, data);
	}
	fclose(fp);
}

static void set_asy_tools(nb_graphics_context_t *g)
{
	g->move_to = nb_graphics_asy_move_to;
	g->line_to = nb_graphics_asy_line_to;
	g->arc = nb_graphics_asy_arc;
	g->set_circle = nb_graphics_asy_set_circle;
	g->set_rectangle = nb_graphics_asy_set_rectangle;
	g->close_path = nb_graphics_asy_close_path;
	g->set_line_width = nb_graphics_asy_line_width;
	g->set_source_rgb = nb_graphics_asy_source_rgb;
	g->set_source_rgba = nb_graphics_asy_source_rgba;
	g->set_source_grad = nb_graphics_asy_source_grad;
	g->set_source_trg = nb_graphics_asy_source_trg;
	g->fill = nb_graphics_asy_fill;
	g->fill_preserve = nb_graphics_asy_fill_preserve;
	g->stroke = nb_graphics_asy_stroke;
	g->set_font_type = nb_graphics_asy_set_font_type;
	g->set_font_size = nb_graphics_asy_set_font_size;
	g->show_text = nb_graphics_asy_show_text;
	g->get_text_attr = nb_graphics_asy_get_text_attr;
}

void nb_graphics_set_camera(nb_graphics_context_t *g, const camera_t *cam)
{
	g->cam = cam;
}

void nb_graphics_unset_camera(nb_graphics_context_t *g)
{
	g->cam = NULL;
}

static void get_camera_view(const camera_t *cam, double cam_vtx[2],
			    double x, double y)
{
	if (NULL != camera) {
		cam_vtx[0] = cam->zoom * (x - cam->center[0]) +
			cam->width / 2.0;
		cam_vtx[1] = cam->zoom * (cam->center[1] - y) +
			cam->height / 2.0;
	} else {
		cam_vtx[0] = x;
		cam_vtx[1] = y;
	}
}

void nb_graphics_move_to(nb_graphics_context_t *g, double x, double y)
{
	double c[2];
	get_camera_view(g->cam, c, x, y);
	g->move_to(g->draw_ptr, c[0], c[1]);
}

void nb_graphics_line_to(nb_graphics_context_t *g, double x, double y)
{
	double c[2];
	get_camera_view(g->cam, c, x, y);
	g->line_to(g->draw_ptr, c[0], c[1]);
}

void nb_graphics_arc(nb_graphics_context_t *g, double x, double y, double r,
		    double a0, double a1)
{
	double c[2];
	get_camera_view(g->cam, c, x, y);
	g->arc(g->draw_ptr, c[0], c[1], r * cam->zoom, a0, a1);
}

void nb_graphics_set_circle(nb_graphics_context_t *g, double x,
			    double y, double r)
{
	double c[2];
	get_camera_view(g->cam, c, x, y);
	g->set_circle(g->draw_ptr, c[0], c[1], r * cam->zoom);
}

void nb_graphics_set_point(nb_graphics_context_t *g,
			   double x, double y, double size)
{
	double c[2];
	get_camera_view(g->cam, c, x, y);
	g->set_circle(g->draw_ptr, c[0], c[1], size);
}

void nb_graphics_set_rectangle(nb_graphics_context_t *g, double x1, double y1,
			       double x2, double y2)
{
	double c1[2], c2[2];
	get_camera_view(g->cam, c1, x1, y1);
	get_camera_view(g->cam, c2, x2, y2);
	g->set_rectangle(g->draw_ptr, c1[0], c1[1], c2[0], c2[1]);
}

void nb_graphics_close_path(nb_graphics_context_t *g)
{
	g->close_path(g->draw_ptr);
}

void nb_graphics_set_line_width(nb_graphics_context_t *g, double w)
{
	g->set_line_width(g->draw_ptr, w);
}

void nb_graphics_set_source_rgb(nb_graphics_context_t *g,
				uint8_t r, uint8_t g, uint8_t b)
{
	g->set_source_rgb(g->draw_ptr, r, g, b);
}

void nb_graphics_set_source_rgba(nb_graphics_context_t *g,
				 uint8_t r, uint8_t g, uint8_t b, uint8_t a)
{
	g->set_source_rgb(g->draw_ptr, r, g, b, a);
}

void nb_graphics_set_source_grad(nb_graphics_context_t *g,
				 nb_graphics_grad_t grad,
				 double x1, double y1,
				 double x2, double y2,
				 nb_palette_t *pat)
{
	g->set_source_grad(g->draw_ptr, grad, x1, y1, x2, y2, pat);
}

void nb_graphics_set_source_trg(nb_graphics_context_t *g,
				double x1, double y1,
				double x2, double y2,
				double x3, double y3,
				uint8_t r1, uint8_t g1, uint8_t b1,
				uint8_t r2, uint8_t g2, uint8_t b2,
				uint8_t r3, uint8_t g3, uint8_t b3)
{
	g->set_source_trg(g->draw_ptr, x1, y1, x2, y2, x3, y3,
			  r1, g1, b1, r2, g2, b2, r3, g3, b3);
}

void nb_graphics_fill(nb_graphics_context_t *g)
{
	g->fill(g->draw_ptr);
}

void nb_graphics_fill_preserve(nb_graphics_context_t *g)
{
	g->fill_preserve(g->draw_ptr);
}

void nb_graphics_stroke(nb_graphics_context_t *g)
{
	g->stroke(g->draw_ptr);
}

void nb_graphics_set_font_type(nb_graphics_context_t *g, const char *type)
{
	g->set_font_type(g->draw_ptr, type);
}

void nb_graphics_set_font_size(nb_graphics_context_t *g, uint16_t size)
{
	g->set_font_size(g->draw_ptr, size);
}

void nb_graphics_show_text(nb_graphics_context_t *g, const char *str)
{
	g->set_show_text(g->draw_ptr, str);
}

void nb_graphics_get_text_attr(nb_graphics_context_t *g, const char *label,
			       nb_text_attr_t *attr)
{
	g->get_text_attr(g->draw_ptr, label, attr);
}
