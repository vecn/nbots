#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <ctype.h>
#include <alloca.h>

#include "nb/math_bot.h"
#include "nb/graphics_bot/drawing_tools.h"

#include "drawing_tools/pix/pix_drawing.h"
#include "drawing_tools/eps/eps_drawing.h"
#include "drawing_tools/asy/asy_drawing.h"

#define MIN(a,b) (((a)<(b))?(a):(b))

#define _SET_RGB(gctx,r,g,b)				\
	do {						\
		(gctx)->set_source_rgb((gctx)->ctx,	\
				       (r), (g), (b));	\
	} while(0)

enum {
	PIX, EPS, ASY, UNKNOWN
};

struct nb_graphics_palette_s {
	/* The palette defines a serie of RGB colors to
	 * colorize values in [0,1]
	 *
	 *  c1    c2       c3         c4  <- RGB colors
	 *   |_____|________|__________|
	 *   0    0.25     0.57        1  <- Tics
	 */
	uint8_t ntics;   /* Number of tics */
	float* tics;   /* Sorted tics in [0,1] */
	uint8_t* rgb;    /* RGB Colors definition */
};

struct nb_graphics_context_s {
	nb_graphics_camera_t cam;
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
	void (*set_source_grad)(void *ctx,
				nb_graphics_grad_t grad,
				double x1, double y1,
				double x2, double y2,
				nb_graphics_palette_t *pal);
	void (*set_source_trg)(void *ctx,
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

static int get_format(const char *filename);
static const char *get_filename_ext(const char *filename);
static void to_lowercase(char *str);
static void full_pipeline(const char* filename, int width, int height,
			  void (*set_tools)(nb_graphics_context_t *g),
			  void (*draw)(nb_graphics_context_t *g, int w, int h,
				       const void *const data),
			  const void *const data);
static void set_pix_tools(nb_graphics_context_t *g);
static void set_eps_tools(nb_graphics_context_t *g);
static void set_asy_tools(nb_graphics_context_t *g);
static void init_gctx(nb_graphics_context_t *g, int width, int height);
static void write_gctx(nb_graphics_context_t *g, const char *filename);
static void clear_gctx(nb_graphics_context_t *g);
static void get_camera_view(const nb_graphics_context_t *g,
			    double cam_vtx[2],
			    double x, double y);

static nb_graphics_palette_t* palette_get_rainbow(void);
static nb_graphics_palette_t* palette_get_sunset(void);
static nb_graphics_palette_t* palette_get_french(void);

static void palette_draw_rectangle(nb_graphics_context_t *g,
				   const nb_graphics_palette_t *const palette,
				   float x, float y, float w, float h,
				   float border);
static void palette_draw_zero_mark(nb_graphics_context_t *g,
				   const nb_graphics_palette_t *const palette,
				   float x, float y, float w, float h,
				   double min_v, double max_v);
static void palette_draw_labels(nb_graphics_context_t *g,
				float font_size, float x, float y,
				float w, float h, double min_v, double max_v);

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
	to_lowercase(ext);
	int format;
	if (0 == strcmp(ext, "png"))
		format = PIX;
	else if (0 == strcmp(ext, "eps"))
		format = EPS;
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

static void to_lowercase(char *str)
{
	int i = 0;
	while (str[i] != '\0' && str[i] != '\n') {
		str[i] = tolower(str[i]);
		i++;
	}
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
	g->set_line_width = nb_graphics_pix_set_line_width;
	g->set_source_rgb = nb_graphics_pix_set_source_rgb;
	g->set_source_rgba = nb_graphics_pix_set_source_rgba;
	g->set_source_grad = nb_graphics_pix_set_source_grad;
	g->set_source_trg = nb_graphics_pix_set_source_trg;
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
	g->set_line_width = nb_graphics_eps_set_line_width;
	g->set_source_rgb = nb_graphics_eps_set_source_rgb;
	g->set_source_rgba = nb_graphics_eps_set_source_rgba;
	g->set_source_grad = nb_graphics_eps_set_source_grad;
	g->set_source_trg = nb_graphics_eps_set_source_trg;
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
	g->set_line_width = nb_graphics_asy_set_line_width;
	g->set_source_rgb = nb_graphics_asy_set_source_rgb;
	g->set_source_rgba = nb_graphics_asy_set_source_rgba;
	g->set_source_grad = nb_graphics_asy_set_source_grad;
	g->set_source_trg = nb_graphics_asy_set_source_trg;
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

nb_graphics_camera_t* nb_graphics_get_camera(nb_graphics_context_t *g)
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
	get_camera_view(g, c, x, y);
	g->move_to(g->ctx, c[0], c[1]);
}

void nb_graphics_line_to(nb_graphics_context_t *g, double x, double y)
{
	double c[2];
	get_camera_view(g, c, x, y);
	g->line_to(g->ctx, c[0], c[1]);
}

void nb_graphics_arc(nb_graphics_context_t *g, double x, double y, double r,
		    double a0, double a1)
{
	double c[2];
	get_camera_view(g, c, x, y);
	g->arc(g->ctx, c[0], c[1], r * g->cam.zoom, a0, a1);
}

void nb_graphics_set_circle(nb_graphics_context_t *g, double x,
			    double y, double r)
{
	double c[2];
	get_camera_view(g, c, x, y);
	g->set_circle(g->ctx, c[0], c[1], r * g->cam.zoom);
}

void nb_graphics_set_point(nb_graphics_context_t *g,
			   double x, double y, double size)
{
	double c[2];
	get_camera_view(g, c, x, y);
	g->set_circle(g->ctx, c[0], c[1], 0.5 * size);
}

void nb_graphics_set_rectangle(nb_graphics_context_t *g, double x1, double y1,
			       double x2, double y2)
{
	double c1[2], c2[2];
	get_camera_view(g, c1, x1, y1);
	get_camera_view(g, c2, x2, y2);
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
		_SET_RGB(g, 0, 0, 0);
		break;
	}
}

void nb_graphics_set_source_rgb(nb_graphics_context_t *gctx,
				uint8_t r, uint8_t g, uint8_t b)
{
	gctx->set_source_rgb(gctx->ctx, r, g, b);
}

void nb_graphics_set_source_rgba(nb_graphics_context_t *gctx,
				 uint8_t r, uint8_t g, uint8_t b, uint8_t a)
{
	gctx->set_source_rgba(gctx->ctx, r, g, b, a);
}

void nb_graphics_set_source_grad(nb_graphics_context_t *g,
				 nb_graphics_grad_t grad,
				 double x1, double y1,
				 double x2, double y2,
				 nb_graphics_palette_t *pal)
{
	g->set_source_grad(g->ctx, grad, x1, y1, x2, y2, pal);
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
	g->show_text(g->ctx, str);
}

void nb_graphics_get_text_attr(nb_graphics_context_t *g, const char *label,
			       nb_text_attr_t *attr)
{
	g->get_text_attr(g->ctx, label, attr);
}
void nb_graphics_cam_fit_box(nb_graphics_camera_t *cam, const double box[4],
			     double width, double height)
{
	cam->center[0] = (box[0] + box[2]) / 2.0;
	cam->center[1] = (box[1] + box[3]) / 2.0;
	cam->zoom = width / (box[2] - box[0]);
	if (cam->zoom > height / (box[3] - box[1]))
		cam->zoom = height / (box[3] - box[1]);
	cam->zoom *= 0.9;
	cam->width = width;
	cam->height = height;
}

nb_graphics_palette_t* nb_graphics_palette_create(void)
{
	return calloc(1, sizeof(nb_graphics_palette_t));
}

nb_graphics_palette_t* nb_graphics_palette_create_preset(nb_graphics_palette_preset preset)
{
	if (NB_RAINBOW == preset)
		return palette_get_rainbow();

	if (NB_SUNSET == preset)
		return palette_get_sunset();

	if (NB_FRENCH == preset)
		return palette_get_french();

	return NULL;
}

void nb_graphics_palette_destroy(nb_graphics_palette_t* palette)
{
	if (palette->ntics > 0) {
		free(palette->tics);
		free(palette->rgb);
	}
	free(palette);
}

void nb_graphics_palette_clear(nb_graphics_palette_t *palette)
{
	if (palette->ntics > 0) {
		free(palette->tics);
		free(palette->rgb);
		palette->tics = NULL;
		palette->rgb = NULL;
	}
	palette->ntics = 0;
}

void nb_graphics_palette_add_colour(nb_graphics_palette_t* palette, float tic,
			    uint8_t r, uint8_t g, uint8_t b)
{
	if (tic < 0)
		tic = 0;
	if (tic > 1)
		tic = 1;
	if (palette->ntics == 0) {
		/* Insert first color */
		palette->tics = (float*)malloc(sizeof(float));
		palette->tics[0] = tic;
		palette->rgb = (uint8_t*)malloc(3);
		palette->rgb[0] = r;
		palette->rgb[1] = g;
		palette->rgb[2] = b;
		palette->ntics = 1;
	} else {
		/* create a new space */
		float* tics = (float*)malloc(palette->ntics*sizeof(float));
		uint8_t* rgb = (uint8_t*)malloc(palette->ntics*3);
		memcpy(tics, palette->tics, palette->ntics*sizeof(float));
		memcpy(rgb, palette->rgb, palette->ntics*3);
		free(palette->tics);
		free(palette->rgb);
		palette->ntics += 1;
		palette->tics = (float*)malloc(palette->ntics*sizeof(float));
		palette->rgb = (uint8_t*)malloc(palette->ntics*3);
		memcpy(palette->tics, tics, (palette->ntics-1)*sizeof(float));
		memcpy(palette->rgb, rgb, (palette->ntics-1)*3);
		free(rgb);
		free(tics);
		/* Insert new color */
		palette->tics[palette->ntics-1] = 2;
		for (uint32_t i=0; i<palette->ntics; i++) {
			if (tic < palette->tics[i]) {
				float aux1 = tic;
				tic = palette->tics[i];
				palette->tics[i] = aux1;
				uint8_t aux2[3] = {r, g, b};
				r = palette->rgb[i * 3];
				g = palette->rgb[i*3+1];
				b = palette->rgb[i*3+2];
				palette->rgb[i * 3] = aux2[0];
				palette->rgb[i*3+1] = aux2[1];
				palette->rgb[i*3+2] = aux2[2];
			}
		}
	}
}

void nb_graphics_palette_get_colour(const nb_graphics_palette_t *const palette,
			    float factor,
			    uint8_t rgb[3])
{
	if (factor <= palette->tics[0]) {
		memcpy(rgb, palette->rgb, 3);
	} else if (factor >= palette->tics[palette->ntics-1]) {
		memcpy(rgb, &(palette->rgb[(palette->ntics-1)*3]), 3);
	} else {
		uint32_t i = 1;
		while(factor > palette->tics[i])
			i++;
		float w1 = (palette->tics[i]-factor)/(palette->tics[i]-palette->tics[i-1]);
		float w2 = (factor-palette->tics[i-1])/(palette->tics[i]-palette->tics[i-1]);
		rgb[0] = w1 * palette->rgb[(i-1) * 3] + w2 * palette->rgb[i * 3];
		rgb[1] = w1 * palette->rgb[(i-1)*3+1] + w2 * palette->rgb[i*3+1];
		rgb[2] = w1 * palette->rgb[(i-1)*3+2] + w2 * palette->rgb[i*3+2];
	}
}

static nb_graphics_palette_t* palette_get_rainbow()
{
	nb_graphics_palette_t* palette = nb_graphics_palette_create();
	nb_graphics_palette_add_colour(palette, 0.00f,   0,   0, 128);
	nb_graphics_palette_add_colour(palette, 0.10f,   0,   0, 255);
	nb_graphics_palette_add_colour(palette, 0.20f,   0, 128, 255);
	nb_graphics_palette_add_colour(palette, 0.37f,   0, 255, 255);
	nb_graphics_palette_add_colour(palette, 0.50f,   0, 255,   0);
	nb_graphics_palette_add_colour(palette, 0.63f, 255, 255,   0);
	nb_graphics_palette_add_colour(palette, 0.80f, 255, 128,   0);
	nb_graphics_palette_add_colour(palette, 0.90f, 255,   0,   0);
	nb_graphics_palette_add_colour(palette, 1.00f, 100,   0,   0);
	return palette;
}

static nb_graphics_palette_t* palette_get_sunset()
{
	nb_graphics_palette_t* palette = nb_graphics_palette_create();
	nb_graphics_palette_add_colour(palette, 0.00f,   0,   0,   0);
	nb_graphics_palette_add_colour(palette, 0.15f,  20,   0, 100);
	nb_graphics_palette_add_colour(palette, 0.30f, 100,   0, 200);
	nb_graphics_palette_add_colour(palette, 0.80f, 220, 100,   0);
	nb_graphics_palette_add_colour(palette, 1.00f, 255, 255,   0);
	return palette;
}

static nb_graphics_palette_t* palette_get_french()
{
	nb_graphics_palette_t* palette = nb_graphics_palette_create();
	nb_graphics_palette_add_colour(palette, 0.00f,   0,   0, 150);
	nb_graphics_palette_add_colour(palette, 0.20f,   0,   0, 255);
	nb_graphics_palette_add_colour(palette, 0.30f, 180, 180, 255);
	nb_graphics_palette_add_colour(palette, 0.50f, 255, 255, 255);
	nb_graphics_palette_add_colour(palette, 0.70f, 255, 180, 180);
	nb_graphics_palette_add_colour(palette, 0.80f, 255,   0,   0);
	nb_graphics_palette_add_colour(palette, 1.00f, 150,   0,   0);
	return palette;
}

void nb_graphics_draw_palette(nb_graphics_context_t *g,
			      const nb_graphics_palette_t *const palette,
			      float x, float y, float w, float h,
			      float border, double min_v, double max_v)
{
	bool cam_status = nb_graphics_is_camera_enabled(g);
	if (cam_status)
		nb_graphics_disable_camera(g);

	palette_draw_rectangle(g, palette, x, y, w, h, border);
	palette_draw_zero_mark(g, palette, x, y, w, h, min_v, max_v);
	palette_draw_labels(g, 10.0f, x, y, w, h, min_v, max_v);

	if (cam_status)
		nb_graphics_enable_camera(g);
}

static void palette_draw_rectangle(nb_graphics_context_t *g,
				   const nb_graphics_palette_t *const palette,
				   float x, float y, float w, float h,
				   float border)
{  
	if (0.0f < border) {
		nb_graphics_set_source(g, NB_BLACK);
		nb_graphics_set_rectangle(g, x-border, y-border,
					  w+2*border, h+2*border);
		nb_graphics_fill(g);
	}

	nb_graphics_set_source_grad(g, NB_LINEAR, x, y+h, x, y, palette);
	nb_graphics_set_rectangle(g, x, y, w, h);
	nb_graphics_fill(g);
}

static void palette_draw_zero_mark(nb_graphics_context_t *g,
				   const nb_graphics_palette_t *const palette,
				   float x, float y, float w, float h,
				   double min_v, double max_v)
{
	if (0 > min_v * max_v) {
		float factor = - min_v / (max_v - min_v);
		uint8_t rgb[3];
		nb_graphics_palette_get_colour(palette, factor, rgb);
	
		rgb[0] = (rgb[0] + 128) % 256;
		rgb[1] = (rgb[1] + 128) % 256;
		rgb[2] = (rgb[2] + 128) % 256;

		nb_graphics_set_source_rgb(g, rgb[0], rgb[1], rgb[2]);

		double yzero = h * factor;
		nb_graphics_move_to(g, x, y + h - yzero);
		nb_graphics_line_to(g, x + w, y + h - yzero);
		nb_graphics_stroke(g);
	}
}

static void palette_draw_labels(nb_graphics_context_t *g,
				float font_size, float x, float y,
				float w, float h, double min_v, double max_v)
{
	nb_graphics_set_font_type(g, "Sans");
	nb_graphics_set_font_size(g, font_size);
	nb_text_attr_t text_attr;
	nb_graphics_set_source(g, NB_BLACK);

	char label[15];
	int n_labels = MIN((int)(h / (font_size + 2)), 10);
	double step_v = (max_v - min_v) / (n_labels - 1.0);
	for (int i = 0; i < n_labels; i++) {
		sprintf(label, "%.3e", max_v - i * step_v);
		nb_graphics_get_text_attr(g, label, &text_attr);

		nb_graphics_move_to(g, x + w + 5.0f, 
				    y + text_attr.height/2.0 + 
				    i * h / (n_labels-1.0));
		nb_graphics_show_text(g, label);
	}
}
