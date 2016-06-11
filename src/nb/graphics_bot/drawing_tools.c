#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <ctype.h>
#include <alloca.h>

#include "nb/math_bot.h"
#include "nb/graphics_bot/drawing_tools.h"

#include "drawing_tools/pix/pix_drawing.h";
#include "drawing_tools/eps/eps_drawing.h";
#include "drawing_tools/asy/asy_drawing.h";

#define _SET_RGB(gctx,r,g,b)				\
	do {						\
		(gctx)->set_source_rgb((gctx)->ctx,	\
				       (r), (g), (b));	\
	} while(0)

enum {
	PIX, EPS, ASY, UNKNOWN;
};

static int get_format(const char *filename);
static const char *get_filename_ext(const char *filename);
static void set_pix_tools(nb_graphics_context_t *g);
static void set_eps_tools(nb_graphics_context_t *g);
static void set_asy_tools(nb_graphics_context_t *g);
static void init_gctx(nb_graphics_context_t *g, int width, int height);
static void write_gctx(nb_graphics_context_t *g, const char *filename);
static void clear_gctx(nb_graphics_context_t *g);
static void get_camera_view(const nb_graphics_context_t *g,
			    double cam_vtx[2],
			    double x, double y);

struct nb_graphics_context_s {
	const nb_graphics_camera_t cam;
	bool using_cam;
	void *ctx;
	void* (*create_context)(int width, int height);
	void (*destroy_context)(void *ctx);
	void (*export_context)(void *ctx, const char *filename);

	void (*move_to)(void *ctx, double x, double y);
	void (*line_to)(void *ctx, double x, double y);
	void (*arc)(void *ctx, double x, double y, double r,
		    double a0, double a1);
	void (*set_circle)(void *ctx, double x, double y, double r);
	void (*set_rectangle)(void *ctx, double x1, double y1,
				  double x2, double y2);
	void (*close_path)(void *ctx);
	void (*set_line_width)(void *ctx, double w);
	void (*set_source_rgb)(void *ctx, uint8_t r, uint8_t g, uint8_t b);
	void (*set_source_rgba)(void *ctx, uint8_t r, uint8_t g,
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
	void (*fill)(void *ctx);
	void (*fill_preserve)(void *ctx);
	void (*stroke)(void *ctx);
	void (*set_font_type)(void *ctx, const char *type);
	void (*set_font_size)(void *ctx, uint16_t size);
	void (*show_text)(void *ctx, const char *str);
	void (*get_text_attr)(void *ctx, const char *label,
			      nb_text_attr_t *attr);
};

void nb_graphics_export(const char* filename, int width, int height,
			void (*draw)(nb_graphics_context_t *g, int w, int h,
				     const void *const data),
			const void *const data)
{
	int format = get_format(filename);
	void *tools;
	switch (format) {
	case PIX:
		tools = set_pix_tools;
		break;
	case EPS:
		tools = set_eps_tools;
		break;
	case ASY:
		tools = set_asy_tools;
		break;
	default:
		tools = set_pix_tools;
	}
	full_pipeline(filename, width, height, tools, draw, data);
}

static int get_format(const char *filename)
{
	const char *raw_ext = get_filename_ext(filename);
	char ext[10];
	strcpy(ext, raw_ext);
	tolower(ext);
	int format;
	if (0 == strcmp(ext, "png"))
		format = PIX;
	else if (0 == strcmp(ext, "eps"))
		format = EPS
	else if (0 == strcmp(ext, "asy"))
		format = ASY;
	else
		format = UNKNOWN;
	return format;
}

static const char *get_filename_ext(const char *filename)
{
	const char *dot = strrchr(filename, '.');
	if (!dot || dot == filename)
		return "";
	else
		return dot + 1;
}

static void full_pipeline(const char* filename, int width, int height,
			  void (*set_tools)(nb_graphics_context_t *g),
			  void (*draw)(nb_graphics_context_t *g, int w, int h,
				       const void *const data),
			  const void *const data)
{
	nb_graphics_context_t g;
	set_tools(&g);
	init_gctx(&g, width, height);

	draw(&g, width, height, data);

	write_gctx(&g, filename);
	clear_gctx(&g);
}

static void set_pix_tools(nb_graphics_context_t *g)
{
	g->create_context = nb_graphics_pix_create_context;
	g->destroy_context = nb_graphics_pix_destroy_context;
	g->export_context = nb_graphics_pix_export_context;
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

static void set_eps_tools(nb_graphics_context_t *g)
{
	g->create_context = nb_graphics_eps_create_context;
	g->destroy_context = nb_graphics_eps_destroy_context;
	g->export_context = nb_graphics_eps_export_context;
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

static void set_asy_tools(nb_graphics_context_t *g)
{
	g->create_context = nb_graphics_asy_create_context;
	g->destroy_context = nb_graphics_asy_destroy_context;
	g->export_context = nb_graphics_asy_export_context;
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

static void init_gctx(nb_graphics_context_t *g, int width, int height)
{
	g->cam.width = width;
	g->cam.height = height;
	g->cam.center[0] = width / 2.0;
	g->cam.center[1] = height / 2.0;
	g->cam.zoom = 1.0;
	g->using_cam = false;
	g->ctx = g->create_context(width, height);
}

static void write_gctx(nb_graphics_context_t *g, const char *filename)
{
	g->export_context(g->ctx, filename);
}

static void clear_gctx(nb_graphics_context_t *g)
{
	g->destroy_context(g->ctx);
	g->ctx = NULL;
}

nb_graphics_camera_t* nb_graphics_get_camera(nb_graphics_context_t *g);
{
	return &(g->cam);
}

bool nb_graphics_is_camera_enabled(const nb_graphics_context_t *g)
{
	return g->using_cam;
}

void nb_graphics_disable_camera(nb_graphics_context_t *g)
{
	g->using_cam = false;
}

void nb_graphics_enable_camera(nb_graphics_context_t *g)
{
	g->using_cam = true;
}

static void get_camera_view(const nb_graphics_context_t *g,
			    double cam_vtx[2],
			    double x, double y)
{
	if (g->using_cam) {
		cam_vtx[0] = g->cam.zoom * (x - g->cam.center[0]) +
			g->cam.width / 2.0;
		cam_vtx[1] = g->cam.zoom * (g->cam.center[1] - y) +
			g->cam.height / 2.0;
	} else {
		cam_vtx[0] = x;
		cam_vtx[1] = y;
	}
}

void nb_graphics_move_to(nb_graphics_context_t *g, double x, double y)
{
	double c[2];
	get_camera_view(g->cam, c, x, y);
	g->move_to(g->ctx, c[0], c[1]);
}

void nb_graphics_line_to(nb_graphics_context_t *g, double x, double y)
{
	double c[2];
	get_camera_view(g->cam, c, x, y);
	g->line_to(g->ctx, c[0], c[1]);
}

void nb_graphics_arc(nb_graphics_context_t *g, double x, double y, double r,
		    double a0, double a1)
{
	double c[2];
	get_camera_view(g->cam, c, x, y);
	g->arc(g->ctx, c[0], c[1], r * cam->zoom, a0, a1);
}

void nb_graphics_set_circle(nb_graphics_context_t *g, double x,
			    double y, double r)
{
	double c[2];
	get_camera_view(g->cam, c, x, y);
	g->set_circle(g->ctx, c[0], c[1], r * cam->zoom);
}

void nb_graphics_set_point(nb_graphics_context_t *g,
			   double x, double y, double size)
{
	double c[2];
	get_camera_view(g->cam, c, x, y);
	g->set_circle(g->ctx, c[0], c[1], 0.5 * size);
}

void nb_graphics_set_rectangle(nb_graphics_context_t *g, double x1, double y1,
			       double x2, double y2)
{
	double c1[2], c2[2];
	get_camera_view(g->cam, c1, x1, y1);
	get_camera_view(g->cam, c2, x2, y2);
	g->set_rectangle(g->ctx, c1[0], c1[1], c2[0], c2[1]);
}

void nb_graphics_close_path(nb_graphics_context_t *g)
{
	g->close_path(g->ctx);
}

void nb_graphics_set_line_width(nb_graphics_context_t *g, double w)
{
	g->set_line_width(g->ctx, w);
}

void nb_graphics_set_source(nb_graphics_context_t *g,
			    nb_graphics_color_t color)
{
	switch(color) {
	case NB_BLACK:
		_SET_RGB(g, 0, 0, 0);
		break;
	case NB_WHITE:
		_SET_RGB(g, 255, 255, 255);
		break;
	case NB_GRAY:
		_SET_RGB(g, 128, 128, 128);
		break;
	case NB_LIGHT_GRAY:
		_SET_RGB(g, 200, 200, 200);
		break;
	case NB_DARK_GRAY:
		_SET_RGB(g, 50, 50, 50);
		break;
	case NB_RED:
		_SET_RGB(g, 255, 0, 0);
		break;
	case NB_GREEN:
		_SET_RGB(g, 0, 255, 0);
		break;
	case NB_BLUE:
		_SET_RGB(g, 0, 0, 255);
		break;
	case NB_YELLOW:
		_SET_RGB(g, 255, 255, 0);
		break;
	case NB_CYAN:
		_SET_RGB(g, 0, 255, 255);
		break;
	case NB_MAGENTA:
		_SET_RGB(g, 255, 0, 255);
		break;
	case NB_ORANGE:
		_SET_RGB(g, 255, 128, 0);
		break;
	case NB_AQUAMARIN:
		_SET_RGB(g, 0, 255, 128);
		break;
	case NB_VIOLET:
		_SET_RGB(g, 128, 0, 255);
		break;
	case NB_ROSE:
		_SET_RGB(g, 255, 0, 128);
		break;
	case NB_CHARTREUSE:
		_SET_RGB(g, 128, 255, 0);
		break;
	case NB_AZURE:
		_SET_RGB(g, 0, 128, 255);
		break;
	default:
		r = 0;
		g = 0;
		b = 0;
		break;
	}
}

void nb_graphics_set_source_rgb(nb_graphics_context_t *g,
				uint8_t r, uint8_t g, uint8_t b)
{
	g->set_source_rgb(g->ctx, r, g, b);
}

void nb_graphics_set_source_rgba(nb_graphics_context_t *g,
				 uint8_t r, uint8_t g, uint8_t b, uint8_t a)
{
	g->set_source_rgb(g->ctx, r, g, b, a);
}

void nb_graphics_set_source_grad(nb_graphics_context_t *g,
				 nb_graphics_grad_t grad,
				 double x1, double y1,
				 double x2, double y2,
				 nb_palette_t *pat)
{
	g->set_source_grad(g->ctx, grad, x1, y1, x2, y2, pat);
}

void nb_graphics_set_source_trg(nb_graphics_context_t *g,
				double x1, double y1,
				double x2, double y2,
				double x3, double y3,
				uint8_t r1, uint8_t g1, uint8_t b1,
				uint8_t r2, uint8_t g2, uint8_t b2,
				uint8_t r3, uint8_t g3, uint8_t b3)
{
	g->set_source_trg(g->ctx, x1, y1, x2, y2, x3, y3,
			  r1, g1, b1, r2, g2, b2, r3, g3, b3);
}

void nb_graphics_fill(nb_graphics_context_t *g)
{
	g->fill(g->ctx);
}

void nb_graphics_fill_preserve(nb_graphics_context_t *g)
{
	g->fill_preserve(g->ctx);
}

void nb_graphics_stroke(nb_graphics_context_t *g)
{
	g->stroke(g->ctx);
}

void nb_graphics_set_font_type(nb_graphics_context_t *g, const char *type)
{
	g->set_font_type(g->ctx, type);
}

void nb_graphics_set_font_size(nb_graphics_context_t *g, uint16_t size)
{
	g->set_font_size(g->ctx, size);
}

void nb_graphics_show_text(nb_graphics_context_t *g, const char *str)
{
	g->set_show_text(g->ctx, str);
}

void nb_graphics_get_text_attr(nb_graphics_context_t *g, const char *label,
			       nb_text_attr_t *attr)
{
	g->get_text_attr(g->ctx, label, attr);
}
