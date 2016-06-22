#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>
#include <alloca.h>
#include <math.h>

#include "nb/image_bot.h"
#include "nb/graphics_bot/drawing_tools.h"

#include "../../palette_struct.h"

#include "rasterizer.h"
#include "font_rasterizer.h"

#include "pix_drawing.h"

#define TURTLE_STATIC_MEMSIZE 20       /* turtle_step */
#define TURTLE_DYNAMIC_MEMINCREASE 25  /* turtle_step */
#define PIXMASK_STATIC_MEMSIZE 2500    /* Bytes */

#define SQRT3 (1.73205080756887729352)

#define ANTIALIASING true

#define POW2(a) ((a)*(a))
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

enum {
	SOLID, GRAD, TRG
};

enum {
	MOVE_TO, LINE_TO, QCURVE_TO,
	QRCURVE_TO, CURVE_TO, CLOSE_PATH
};

typedef struct {
	uint8_t type;
	float x, y;
	float data[4];
} turtle_step;

typedef struct {
	uint16_t N;
	uint16_t N_alloc;
	turtle_step static_mem[TURTLE_STATIC_MEMSIZE];
	turtle_step *dynamic_mem;
} turtle_t;

typedef struct {
	uint8_t source_type:4;
	nb_graphics_grad_t grad_type:4;
	float vtx[6];
	nb_graphics_palette_t *pal;
	nb_graphics_palette_t color;
} source_t;

typedef struct {
	const char *type;
	uint16_t size;	
} font_t;

typedef struct {
	int32_t xmin, ymin;
	int32_t width, height;
	uint8_t static_pix[PIXMASK_STATIC_MEMSIZE];
	uint8_t *pix;
} rasterized_path_t;

typedef struct {
	vcn_image_t *img;
	turtle_t *turtle;
	rasterized_path_t *pen_stencil;
	float line_width;
	source_t *source;
	font_t *font;
} context_t;

static void turtle_clear(turtle_t *turtle);
static void turtle_add(turtle_t *turtle, uint8_t type,
		       float x, float y, float v1, float v2,
		       float v3, float v4);
static void set_pen_stencil(rasterized_path_t *pen_stencil, float thickness);
static void set_pen_stencil_w2(rasterized_path_t *pen_stencil, float thickness);
static void set_pen_stencil_w3(rasterized_path_t *pen_stencil, float thickness);
static void set_pen_stencil_w4(rasterized_path_t *pen_stencil, float thickness);
static void set_turtle_stencil(turtle_t *turtle_stencil, float thickness);
static void set_line_pixel(int x, int y, uint8_t i, void* context);
static void set_pen_pixel(int x, int y, uint8_t i, void* context);
static turtle_step *turtle_ref_step(turtle_t *turtle, uint16_t i);
static void turtle_reset(turtle_t *turtle);
static void set_pixel(int x, int y, uint8_t i, void* context);
static void draw_line(context_t *c, int x0, int y0, int x1, int y1);
static void draw_qcurve(context_t *c, int x0, int y0, int x1, int y1,
			int cx, int cy);
static void draw_qrcurve(context_t *c, int x0, int y0, int x1, int y1,
			 int cx, int cy, float w);
static void draw_curve(context_t *c, int x0, int y0, int x1, int y1,
		       float c0x, float c0y, float c1x, float c1y);
static void source_set_rgba(source_t *source, uint8_t r,
			    uint8_t g, uint8_t b, uint8_t a);
static void source_set_grad(source_t *source,
			    nb_graphics_grad_t grad,
			    float x1, float y1,
			    float x2, float y2,
			    nb_graphics_palette_t *pal);
static void source_set_trg(source_t *source,
			   float x1, float y1,
			   float x2, float y2,
			   float x3, float y3,
			   const uint8_t rgba1[4],
			   const uint8_t rgba2[4],
			   const uint8_t rgba3[4]);
static void rpath_init(rasterized_path_t *rpath, const turtle_t *turtle);
static void rpath_alloc_pix(rasterized_path_t *rpath);
static void turtle_step_get_box(const turtle_step *step, uint32_t box[4]);
static void rpath_full(rasterized_path_t *rpath);
static void rpath_unfill_outside(rasterized_path_t *rpath);
static void rpath_set_pixel(int x, int y, uint8_t i, void *rpath);
static void rpath_unset_pixel(int x, int y, uint8_t i, void *rpath);
static uint8_t rpath_get_intensity(const void *rpath, int x, int y);
static bool rpath_pixel_is_not_empty(int x, int y, const void *rpath);
static void rpath_blend_image(const rasterized_path_t *rpath,
			      const source_t *source,
			      uint8_t intensity,
			      vcn_image_t *img);
static void rpath_finish(rasterized_path_t *rpath);
static void rasterize_turtle(const turtle_t *turtle,
			     bool antialiased,
			     void (*set_pixel)(int x, int y, uint8_t i, void *),
			     void *pixel_data);
static void source_get_color(const source_t *source, int x, int y,
			     uint8_t pix[4]);
static void source_get_color_solid(const source_t *source, uint8_t pix[4]);
static void source_get_color_grad(const source_t *source, int x, int y,
				  uint8_t pix[4]);
static void source_get_color_grad_linear(const source_t *source, int x, int y,
					 uint8_t pix[4]);
static void source_get_color_grad_radial(const source_t *source, int x, int y,
					 uint8_t pix[4]);
static void source_get_color_trg(const source_t *source, int x, int y,
				 uint8_t pix[4]);
static void get_barycentric_coordinates(float x1, float y1, float x2, float y2,
					float x3, float y3, float x, float y,
					float lambda[3]);

void* nb_graphics_pix_create_context(int width, int height)
{
	uint32_t ctx_size = sizeof(context_t);
	uint32_t img_size = sizeof(vcn_image_t);
	uint32_t pix_size = 4 * width * height;
	uint32_t trt_size = sizeof(turtle_t);
	uint16_t rp_size = sizeof(rasterized_path_t);
	uint16_t src_size = sizeof(source_t);
	uint16_t fnt_size = sizeof(font_t);
	uint32_t memsize = ctx_size + img_size + pix_size +
		trt_size + rp_size + src_size + fnt_size;
	char *memblock = malloc(memsize);
	context_t *ctx = (void*) memblock;

	ctx->img = (void*) (memblock + ctx_size);
	vcn_image_init(ctx->img);
	ctx->img->width = width;
	ctx->img->height = height;
	ctx->img->comp_x_pixel = 4;
	ctx->img->pixels = (void*) (memblock + ctx_size + img_size);
	memset(ctx->img->pixels, 0, pix_size);

	ctx->turtle = (void*) (memblock + ctx_size + img_size + pix_size);
	memset(ctx->turtle, 0, trt_size);

	ctx->pen_stencil = (void*) (memblock + ctx_size + img_size +
				    pix_size + trt_size);
	memset(ctx->pen_stencil, 0, rp_size);

	ctx->line_width = 1.0;
	
	ctx->source = (void*) (memblock + ctx_size + img_size +
			       pix_size + trt_size + rp_size);
	source_set_rgba(ctx->source, 0, 0, 0, 255);
	
	ctx->font = (void*) (memblock + ctx_size + img_size +
			     pix_size + trt_size + rp_size + src_size);
	ctx->font->type = "Mono";
	ctx->font->size = 18;

	return ctx;
}

void nb_graphics_pix_destroy_context(void *ctx)
{
	context_t *c = ctx;
	turtle_clear(c->turtle);
	rpath_finish(c->pen_stencil);
	free(ctx);
}

static void turtle_clear(turtle_t *turtle)
{
	if (turtle->N_alloc > 0)
		free(turtle->dynamic_mem);
	memset(turtle, 0, sizeof(turtle_t));
}

void nb_graphics_pix_export_context(const void *ctx, const char *filename)
{
	const context_t *c = ctx;
	vcn_image_write(c->img, filename);
}

void nb_graphics_pix_move_to(void *ctx, float x, float y)
{
	context_t *c = ctx;
	turtle_add(c->turtle, MOVE_TO, x, y,
		   0, 0, 0, 0);
}

static void turtle_add(turtle_t *turtle, uint8_t type,
		       float x, float y, float v1, float v2,
		       float v3, float v4)
{
	if (turtle->N < TURTLE_STATIC_MEMSIZE) {
		turtle->static_mem[turtle->N].type = type;
		turtle->static_mem[turtle->N].x = x;
		turtle->static_mem[turtle->N].y = y;
		turtle->static_mem[turtle->N].data[0] = v1;
		turtle->static_mem[turtle->N].data[1] = v2;
		turtle->static_mem[turtle->N].data[2] = v3;
		turtle->static_mem[turtle->N].data[3] = v4;
		turtle->N += 1;
	} else {
		uint16_t N = turtle->N - TURTLE_STATIC_MEMSIZE;
		if (N >= turtle->N_alloc) {
			turtle->N_alloc += TURTLE_DYNAMIC_MEMINCREASE;
			uint16_t memsize =
				sizeof(turtle_step) * turtle->N_alloc;
			turtle->dynamic_mem = realloc(turtle->dynamic_mem,
						      memsize);
		}
		turtle->dynamic_mem[N].type = type;
		turtle->dynamic_mem[N].x = x;
		turtle->dynamic_mem[N].y = y;
		turtle->dynamic_mem[N].data[0] = v1;
		turtle->dynamic_mem[N].data[1] = v2;
		turtle->dynamic_mem[N].data[2] = v3;
		turtle->dynamic_mem[N].data[3] = v4;
		turtle->N += 1;
	}
}

void nb_graphics_pix_line_to(void *ctx, float x, float y)
{
	context_t *c = ctx;
	turtle_add(c->turtle, LINE_TO, x, y,
		   0, 0, 0, 0);
}

void nb_graphics_pix_qcurve_to(void *ctx, float x, float y,
			       float xcontrol, float ycontrol)
{
	context_t *c = ctx;
	turtle_add(c->turtle, QCURVE_TO, x, y,
		   xcontrol, ycontrol, 0, 0);
}

void nb_graphics_pix_qrcurve_to(void *ctx, float x, float y,
				float xcontrol, float ycontrol, float w)
{
	context_t *c = ctx;
	turtle_add(c->turtle, QRCURVE_TO, x, y,
		   xcontrol, ycontrol, w, 0);
}
void nb_graphics_pix_curve_to(void *ctx, float x, float y,
			      float x0_control, float y0_control,
			      float x1_control, float y1_control)
{
	context_t *c = ctx;
	turtle_add(c->turtle, CURVE_TO, x, y,
		   x0_control, y0_control,
		   x1_control, y1_control);
}

void nb_graphics_pix_close_path(void *ctx)
{
	context_t *c = ctx;
	turtle_add(c->turtle, CLOSE_PATH,
		   0, 0, 0, 0, 0, 0);
}

void nb_graphics_pix_set_line_width(void *ctx, float w)
{
	context_t *c = ctx;
	c->line_width = fabs(w);
	if (c->line_width > 1.0)
		set_pen_stencil(c->pen_stencil, c->line_width);
}

static void set_pen_stencil(rasterized_path_t *pen_stencil, float thickness)
{
	int w = (int)thickness;
	switch (w) {
	case 1:
		set_pen_stencil_w2(pen_stencil, thickness);
		break;
	case 2:
		set_pen_stencil_w3(pen_stencil, thickness);
		break;
	default:
		set_pen_stencil_w4(pen_stencil, thickness);
	}
}

static void set_pen_stencil_w2(rasterized_path_t *pen_stencil, float thickness)
{
	thickness -= 1.0;
	pen_stencil->width = 2;
	pen_stencil->height = 2;
	pen_stencil->pix = pen_stencil->static_pix;
	pen_stencil->pix[0] = (int)(100.0f * thickness + 0.5);
	pen_stencil->pix[1] = (int)(255.0f * thickness + 0.5);
	pen_stencil->pix[2] = (int)(255.0f * thickness + 0.5);
	pen_stencil->pix[3] = 255;
}

static void set_pen_stencil_w3(rasterized_path_t *pen_stencil, float thickness)
{
	thickness -= 2.0;
	pen_stencil->width = 3;
	pen_stencil->height = 3;
	pen_stencil->pix = pen_stencil->static_pix;
	pen_stencil->pix[0] = (int)(50.0f * thickness + 0.5);
	pen_stencil->pix[1] = (int)(255.0f * thickness + 0.5);
	pen_stencil->pix[2] = (int)(100.0f * thickness + 0.5);
	pen_stencil->pix[3] = (int)(255.0f * thickness + 0.5);
	pen_stencil->pix[4] = 255;
	pen_stencil->pix[5] = 255;
	pen_stencil->pix[6] = (int)(100.0f * thickness + 0.5);
	pen_stencil->pix[7] = 255;
	pen_stencil->pix[8] = 255;
}

static void set_pen_stencil_w4(rasterized_path_t *pen_stencil, float thickness)
{	
	turtle_t *turtle_stencil = alloca(sizeof(turtle_t));
	memset(turtle_stencil, 0, sizeof(turtle_t));
	set_turtle_stencil(turtle_stencil, thickness);
	
	pen_stencil->width = (int)thickness;
	pen_stencil->height = (int)thickness;
	rpath_alloc_pix(pen_stencil);

	rpath_full(pen_stencil);
	rasterize_turtle(turtle_stencil, false, rpath_unset_pixel, pen_stencil);
	rpath_unfill_outside(pen_stencil);
	rasterize_turtle(turtle_stencil, false, rpath_set_pixel, pen_stencil);

	turtle_clear(turtle_stencil);
}

static void set_turtle_stencil(turtle_t *turtle_stencil, float thickness)
{
	float r = (int)((thickness-1) / 2.0f);
	turtle_add(turtle_stencil, MOVE_TO, r, 0.0f, 0, 0, 0, 0);
	turtle_add(turtle_stencil, QRCURVE_TO, r + r * SQRT3/2.0f, r + r/2.0f,
		   r + r * SQRT3, 0.0f, 0.5f, 0);
	turtle_add(turtle_stencil, QRCURVE_TO, r - r * SQRT3/2.0f, r + r/2.0f,
		   r, r + 2.0f * r, 0.5f, 0);
	turtle_add(turtle_stencil, QRCURVE_TO, r, 0.0f,
		   r - r * SQRT3, 0.0f, 0.5f, 0);
}

void nb_graphics_pix_set_source_rgb(void *ctx,
				    uint8_t r, uint8_t g, uint8_t b)
{
	context_t *c = ctx;
	source_set_rgba(c->source, r, g, b, 255);
}

void nb_graphics_pix_set_source_rgba(void *ctx, uint8_t r, uint8_t g,
				     uint8_t b, uint8_t a)
{
	context_t *c = ctx;
	source_set_rgba(c->source, r, g, b, a);
}

static void source_set_rgba(source_t *source, uint8_t r,
			    uint8_t g, uint8_t b, uint8_t a)
{
	memset(source, 0, sizeof(source_t));
	source->source_type = SOLID;
	source->pal = &(source->color);
	nb_graphics_palette_add_rgba(source->pal, 0.0,
				     r, g, b, a);
}
void nb_graphics_pix_set_source_grad(void *ctx,
				     nb_graphics_grad_t grad,
				     float x1, float y1,
				     float x2, float y2,
				     nb_graphics_palette_t *pal)
{
	context_t *c = ctx;
	source_set_grad(c->source, grad, x1, y1, x2, y2, pal);
}

static void source_set_grad(source_t *source,
			    nb_graphics_grad_t grad,
			    float x1, float y1,
			    float x2, float y2,
			    nb_graphics_palette_t *pal)
{
	memset(source, 0, sizeof(source_t));
	source->source_type = GRAD;
	source->grad_type = grad;
	source->pal = pal;
	source->vtx[0] = x1;
	source->vtx[1] = y1;
	source->vtx[2] = x2 - x1;
	source->vtx[3] = y2 - y1;
	source->vtx[4] = sqrt(POW2(source->vtx[2]) + 
			      POW2(source->vtx[3]));
	source->vtx[2] /= source->vtx[4];
	source->vtx[3] /= source->vtx[4];
}

void nb_graphics_pix_set_source_trg(void *ctx,
				    float x1, float y1,
				    float x2, float y2,
				    float x3, float y3,
				    const uint8_t rgba1[4],
				    const uint8_t rgba2[4],
				    const uint8_t rgba3[4])
{
	context_t *c = ctx;
	source_set_trg(c->source, x1, y1, x2, y2, x3, y3,
		       rgba1, rgba2, rgba3);
}

static void source_set_trg(source_t *source,
			   float x1, float y1,
			   float x2, float y2,
			   float x3, float y3,
			   const uint8_t rgba1[4],
			   const uint8_t rgba2[4],
			   const uint8_t rgba3[4])
{
	memset(source, 0, sizeof(source_t));
	source->source_type = TRG;
	source->pal = &(source->color);
	source->vtx[0] = x1;
	source->vtx[1] = y1;
	source->vtx[2] = x2;
	source->vtx[3] = y2;
	source->vtx[4] = x3;
	source->vtx[5] = y3;
	nb_graphics_palette_add_rgba(source->pal, 0.0, rgba1[0],
				     rgba1[1], rgba1[2], rgba1[3]);
	nb_graphics_palette_add_rgba(source->pal, 0.5, rgba2[0],
				     rgba2[1], rgba2[2], rgba2[3]);
	nb_graphics_palette_add_rgba(source->pal, 1.0, rgba3[0],
				     rgba3[1], rgba3[2], rgba3[3]);
}

void nb_graphics_pix_fill(void *ctx)
{
	context_t *c = ctx;
	nb_graphics_pix_fill_preserve(ctx);
	turtle_reset(c->turtle);
}

void nb_graphics_pix_fill_preserve(void *ctx)
{
	
	context_t *c = ctx;
	rasterized_path_t *rpath = alloca(sizeof(rasterized_path_t));
	rpath_init(rpath, c->turtle);
	if (rpath->xmin < c->img->width &&
	    rpath->ymin < c->img->height &&
	    rpath->xmin + rpath->width >= 0 &&
	    rpath->ymin + rpath->height >= 0) {
		rpath_full(rpath);
		rasterize_turtle(c->turtle, false, rpath_unset_pixel, rpath);
		rpath_unfill_outside(rpath);
		rasterize_turtle(c->turtle, true, rpath_set_pixel, rpath);
		rpath_blend_image(rpath, c->source, 255, c->img);
	}
	rpath_finish(rpath);
}

static void rpath_init(rasterized_path_t *rpath, const turtle_t *turtle)
{
	rpath->xmin = 0;
	rpath->ymin = 0;
	uint32_t xmax = 0;
	uint32_t ymax = 0;
	for (uint16_t i = 0; i < turtle->N; i++) {
		turtle_step *step = turtle_ref_step((turtle_t*)turtle, i);
		uint32_t box[4];
		turtle_step_get_box(step, box);

		if (0 == i) {
			rpath->xmin = box[0];
			rpath->ymin = box[1];
			xmax = box[2];
			ymax = box[3];
		} else {
			if (box[0] < rpath->xmin)
				rpath->xmin = box[0];
			else if (box[2] > xmax)
				xmax = box[2];

			if (box[1] < rpath->ymin)
				rpath->ymin = box[1];
			else if (box[3] > ymax)
				ymax = box[3];
		}
	}
	rpath->width = xmax - rpath->xmin + 1;
	rpath->height = ymax - rpath->ymin + 1;
	
	rpath_alloc_pix(rpath);
}

static void rpath_alloc_pix(rasterized_path_t *rpath)
{
	if (rpath->width * rpath->height > PIXMASK_STATIC_MEMSIZE) {
		uint32_t memsize = (uint32_t)(rpath->width * rpath->height);
		rpath->pix = calloc(memsize, 1);
	} else {
		memset(rpath->static_pix, 0, PIXMASK_STATIC_MEMSIZE);
		rpath->pix = rpath->static_pix;
	}
}

static void turtle_step_get_box(const turtle_step *step, uint32_t box[4])
{
	int sx = (int)(step->x + 0.5);
	int sy = (int)(step->y + 0.5);
	if (QCURVE_TO == step->type ||
	    QRCURVE_TO == step->type) {
		int cx = (int)(step->data[0] + 0.5);
		int cy = (int)(step->data[1] + 0.5);
			box[0] = MIN(sx, cx);
			box[2] = MAX(sx, cx);
			box[1] = MIN(sy, cy);
			box[3] = MAX(sy, cy);
	} else if (CURVE_TO) {
		int c0x = (int)(step->data[0] + 0.5);
		int c0y = (int)(step->data[1] + 0.5);
		int c1x = (int)(step->data[0] + 0.5);
		int c1y = (int)(step->data[1] + 0.5);
		if (c0x < c1x) {
			box[0] = MIN(sx, c0x);
			box[2] = MAX(sx, c1x);
		} else {
			box[0] = MIN(sx, c1x);
			box[2] = MAX(sx, c0x);
		}
		if (c0y < c1y) {
			box[1] = MIN(sy, c0y);
			box[3] = MAX(sy, c1y);
		} else {
			box[1] = MIN(sy, c1y);
			box[3] = MAX(sy, c0y);
		}
	} else {
		box[0] = sx;
		box[1] = sy;
		box[2] = sx;
		box[3] = sy;
	}
}

static void rpath_full(rasterized_path_t *rpath)
{
	memset(rpath->pix, -1, rpath->width * rpath->height);
}

static void rpath_unfill_outside(rasterized_path_t *rpath)
{
	for (uint32_t i = 0; i < rpath->width; i++) {
		uint32_t x = rpath->xmin + i;
		uint32_t y =  0;
		nb_graphics_rasterizer_fill(x, y, 0,
					    rpath_unset_pixel,
					    rpath_pixel_is_not_empty,
					    rpath->width,
					    rpath->height,
					    rpath);
		
		y = rpath->height - 1;
		nb_graphics_rasterizer_fill(x, y, 0,
					    rpath_unset_pixel,
					    rpath_pixel_is_not_empty,
					    rpath->width,
					    rpath->height,
					    rpath);
	}
	for (uint32_t i = 1; i < rpath->height - 1; i++) {
		uint32_t x =  0;
		uint32_t y = rpath->ymin + i;
		nb_graphics_rasterizer_fill(x, y, 0,
					    rpath_unset_pixel,
					    rpath_pixel_is_not_empty,
					    rpath->width,
					    rpath->height,
					    rpath);
		x = rpath->width - 1;		
		nb_graphics_rasterizer_fill(x, y, 0,
					    rpath_unset_pixel,
					    rpath_pixel_is_not_empty,
					    rpath->width,
					    rpath->height,
					    rpath);
	}
}

static void rpath_set_pixel(int x, int y, uint8_t i, void *rpath)
{
	rasterized_path_t *rp = rpath;
	x = x - rp->xmin;
	y = y - rp->ymin;
	if (x >= 0 && x < rp->width && y >= 0 && y < rp->height) {
		uint32_t p = y * rp->width + x;
		rp->pix[p] = MAX(rp->pix[p], i);
	}
}

static void rpath_unset_pixel(int x, int y, uint8_t i, void *rpath)
{
	rasterized_path_t *rp = rpath;
	x = x - rp->xmin;
	y = y - rp->ymin;
	if (x >= 0 && x < rp->width && y >= 0 && y < rp->height) {
		uint32_t p = y * rp->width + x;
		rp->pix[p] = 0;
	}
}

static uint8_t rpath_get_intensity(const void *rpath, int x, int y)
{
	const rasterized_path_t *rp = rpath;
	uint8_t out = 0;
	x = x - rp->xmin;
	y = y - rp->ymin;
	if (x >= 0 && x < rp->width && y >= 0 && y < rp->height) {
		uint32_t p = y * rp->width + x;
		out = rp->pix[p];
	}
	return out;
}

static bool rpath_pixel_is_not_empty(int x, int y, const void *rpath)
{
	const rasterized_path_t *rp = rpath;
	bool out = false;
	x = x - rp->xmin;
	y = y - rp->ymin;
	if (x >= 0 && x < rp->width && y >= 0 && y < rp->height) {
		uint32_t p = y * rp->width + x;
		out = (0 < rp->pix[p]);
	}
	return out;
}

static void rpath_blend_image(const rasterized_path_t *rpath,
			      const source_t *source,
			      uint8_t intensity,
			      vcn_image_t *img)
{
	uint8_t pix[4];
	int i0 = MAX(0, -rpath->xmin);
	int j0 = MAX(0, -rpath->ymin);
	int w = MIN(rpath->width, img->width - rpath->xmin);
	int h = MIN(rpath->height, img->height - rpath->ymin);
	for (uint32_t i = i0; i < w; i++) {
		for (uint32_t j = j0; j < h; j++) {
			int32_t x = rpath->xmin + i;
			int32_t y = rpath->ymin + j;
			if (rpath_pixel_is_not_empty(x, y, rpath)) {
				float I = (intensity/255.0f) *
					(rpath_get_intensity(rpath, x, y)/
					 255.0f);
				source_get_color(source, x, y, pix);
				if (I < 1.0)
					pix[3] = (int)(pix[3] * I + 0.5);
				vcn_image_blend_pixel_rgba(img, y, x, pix);
			}
		}
	}
}

static void rpath_finish(rasterized_path_t *rpath)
{
	if (rpath->pix != rpath->static_pix &&
	    rpath->pix != NULL)
		free(rpath->pix);
}

static void rasterize_turtle(const turtle_t *turtle,
			     bool antialiased,
			     void (*set_pixel)(int x, int y, uint8_t i, void *),
			     void *pixel_data)
{
	int x0, y0;
	int x, y;
	for (uint16_t i = 0; i < turtle->N; i++) {
		turtle_step *step = turtle_ref_step(turtle, i);
		int sx = (int)(step->x + 0.5);
		int sy = (int)(step->y + 0.5);
		switch (step->type) {
		case MOVE_TO:
			x0 = sx;
			y0 = sy;
			break;
		case LINE_TO:
			nb_graphics_rasterizer_line(x, y, sx, sy,
						    antialiased,
						    set_pixel,
						    pixel_data);
			break;
		case QCURVE_TO:
			nb_graphics_rasterizer_quad_bezier(x, y, sx, sy,
							   step->data[0] + 0.5,
							   step->data[1] + 0.5,
							   antialiased,
							   set_pixel,
							   pixel_data);
			break;
		case QRCURVE_TO:
			nb_graphics_rasterizer_quad_rational_bezier
				(x, y, sx, sy,
				 step->data[0] + 0.5,
				 step->data[1] + 0.5,
				 step->data[2],
				 antialiased,
				 set_pixel,
				 pixel_data);
			break;
		case CURVE_TO:
			nb_graphics_rasterizer_cubic_bezier(x, y, sx, sy,
							    step->data[0],
							    step->data[1],
							    step->data[2],
							    step->data[3],
							    antialiased,
							    set_pixel,
							    pixel_data);
			break;
		case CLOSE_PATH:
			sx = x0;
			sy = y0;
			nb_graphics_rasterizer_line(x, y, sx, sy,
						    antialiased,
						    set_pixel,
						    pixel_data);
			break;
		default:
			sx = x0;
			sy = y0;
			nb_graphics_rasterizer_line(x, y, sx, sy,
						    antialiased,
						    set_pixel,
						    pixel_data);
		}
		x = sx;
		y = sy;
	}
}

void nb_graphics_pix_stroke(void *ctx)
{
	context_t *c = ctx;
	nb_graphics_pix_stroke_preserve(ctx);
	turtle_reset(c->turtle);
}

void nb_graphics_pix_stroke_preserve(void *ctx)
{
	context_t *c = ctx;
	if (c->line_width <= 1.0)
		rasterize_turtle(c->turtle, true, set_line_pixel, c);
	else
		rasterize_turtle(c->turtle, true, set_pen_pixel, c);
}

static void set_line_pixel(int x, int y, uint8_t i, void* context)
{
	context_t *c = context;
	float w = c->line_width;
	if (x >= 0 && x < c->img->width && y >= 0 && y < c->img->height) {
		uint8_t pix[4];
		source_get_color(c->source, x, y, pix);
		float I = w * (pix[3]/255.0f) * (i/255.0f);
		pix[3] = (int)(255.0f * I + 0.5);
		if (pix[3] > 0)
			vcn_image_blend_pixel_rgba(c->img, y, x, pix);
	}
}

static void set_pen_pixel(int x, int y, uint8_t i, void* context)
{
	context_t *c = context;
	int w_div_2 = (int)(c->pen_stencil->width / 2.0);
	int h_div_2 = (int)(c->pen_stencil->height / 2.0);
	if (x + w_div_2 >= 0 && x - w_div_2 < c->img->width &&
	    y + h_div_2 >= 0 && y - h_div_2 < c->img->height) {
		c->pen_stencil->xmin = x - w_div_2;
		c->pen_stencil->ymin = y - h_div_2;
		rpath_blend_image(c->pen_stencil, c->source, i, c->img);
		c->pen_stencil->xmin = 0;
		c->pen_stencil->ymin = 0;
	}
}

static turtle_step *turtle_ref_step(turtle_t *turtle, uint16_t i)
{
	void *step;
	if (i < TURTLE_STATIC_MEMSIZE)
		step = &(turtle->static_mem[i]);
	else
		step = &(turtle->dynamic_mem[i - TURTLE_STATIC_MEMSIZE]);
	return step;
}

static void turtle_reset(turtle_t *turtle)
{
	turtle->N = 0;
}

static void set_pixel(int x, int y, uint8_t i, void* context)
{
	context_t *c = context;
	if (x >= 0 && x < c->img->width && y >= 0 && y < c->img->height) {
		uint8_t pix[4];
		source_get_color(c->source, x, y, pix);
		pix[3] = (uint8_t)(255.0f * (pix[3]/255.0f) * (i/255.0f) + 0.5);
		if (pix[3] > 0)
			vcn_image_blend_pixel_rgba(c->img, y, x, pix);
	}
}

static void source_get_color(const source_t *source, int x, int y, 
			     uint8_t pix[4])
{
	switch (source->source_type) {
	case SOLID:
		source_get_color_solid(source, pix);
		break;
	case GRAD:
		source_get_color_grad(source, x, y, pix);
		break;
	case TRG:
		source_get_color_trg(source, x, y, pix);
		break;
	default:
		source_get_color_solid(source, pix);
	}
}

static void source_get_color_solid(const source_t *source, uint8_t pix[4])
{
	nb_graphics_palette_get_rgba(source->pal, 0.0, pix);
}

static void source_get_color_grad(const source_t *source, int x, int y,
				  uint8_t pix[4])
{
	switch (source->grad_type) {
	case NB_LINEAR:
		source_get_color_grad_linear(source, x, y, pix);
		break;
	case NB_RADIAL:
		source_get_color_grad_radial(source, x, y, pix);
		break;
	default:
		source_get_color_grad_linear(source, x, y, pix);
	}
}

static void source_get_color_grad_linear(const source_t *source, int x, int y,
					 uint8_t pix[4])
{
	x = x - source->vtx[0];
	y = y - source->vtx[1];
	float factor = x * source->vtx[2] + y * source->vtx[3];
	factor = factor / source->vtx[4];
	nb_graphics_palette_get_rgba(source->pal, factor, pix);
}

static void source_get_color_grad_radial(const source_t *source, int x, int y,
					 uint8_t pix[4])
{
	x = x - source->vtx[0];
	y = y - source->vtx[1];
	float factor = sqrt(POW2(x) + POW2(y));
	factor = factor / source->vtx[4];
	nb_graphics_palette_get_rgba(source->pal, factor, pix);
}

static void source_get_color_trg(const source_t *source, int x, int y,
				 uint8_t pix[4])
{
	float lambda[3];
	get_barycentric_coordinates(source->vtx[0], source->vtx[1],
				    source->vtx[2], source->vtx[3],
				    source->vtx[4], source->vtx[5],
				    x, y, lambda);
	uint8_t pix1[4], pix2[4], pix3[4];
	nb_graphics_palette_get_rgba(source->pal, 0.0, pix1);
	nb_graphics_palette_get_rgba(source->pal, 0.5, pix2);
	nb_graphics_palette_get_rgba(source->pal, 1.0, pix3);

	for (int j = 0; j < 4; j++)
		pix[j] = (uint8_t) (lambda[0] * pix1[j] +
				    lambda[1] * pix2[j] +
				    lambda[2] * pix3[j]);
}

static void get_barycentric_coordinates(float x1, float y1, float x2, float y2,
					float x3, float y3, float x, float y,
					float lambda[3])
{
	float detT = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3);
	float x_x3 = x - x3;
	float y_y3 = y - y3;
	lambda[0] = ((y2 - y3) * x_x3 + (x3 - x2) * y_y3) / detT;
	lambda[1] = ((y3 - y1) * x_x3 + (x1 - x3) * y_y3) / detT;
	lambda[2] = 1.0f - lambda[1] - lambda[0];
	if (lambda[0] < 0) {
			lambda[0] = 0;
			float d2 = sqrt(POW2(x - x2) + POW2(y - y2));
			float d3 = sqrt(POW2(x - x3) + POW2(y - y3));
			lambda[1] = d3 / (d2 + d3);
			lambda[2] = 1.0f - lambda[1];
	} else if (lambda[1] < 0) {
			lambda[1] = 0;
			float d2 = sqrt(POW2(x - x1) + POW2(y - y1));
			float d3 = sqrt(POW2(x - x3) + POW2(y - y3));
			lambda[0] = d3 / (d2 + d3);
			lambda[2] = 1.0f - lambda[0];
	} else if (lambda[2] < 0) {
			lambda[2] = 0;
			float d2 = sqrt(POW2(x - x1) + POW2(y - y1));
			float d3 = sqrt(POW2(x - x2) + POW2(y - y2));
			lambda[0] = d3 / (d2 + d3);
			lambda[1] = 1.0f - lambda[0];
	} 
}

static void draw_line(context_t *c, int x0, int y0, int x1, int y1)
{
	nb_graphics_rasterizer_line(x0, y0, x1, y1,
				    ANTIALIASING,
				    set_pixel, c);
}

static void draw_qcurve(context_t *c, int x0, int y0, int x1, int y1,
			int cx,	int cy)
{
	nb_graphics_rasterizer_quad_bezier(x0, y0, x1, y1, cx, cy,
					   ANTIALIASING,
					   set_pixel, c);
}

static void draw_qrcurve(context_t *c, int x0, int y0, int x1, int y1,
			 int cx, int cy, float w)
{
	nb_graphics_rasterizer_quad_rational_bezier(x0, y0, x1, y1,
						    cx, cy, w,
						    ANTIALIASING,
						    set_pixel, c);
}

static void draw_curve(context_t *c, int x0, int y0, int x1, int y1,
		       float c0x, float c0y, float c1x, float c1y)
{
	nb_graphics_rasterizer_cubic_bezier(x0, y0, x1, y1,
					    c0x, c0y, c1x, c1y,
					    ANTIALIASING,
					    set_pixel, c);
}

void nb_graphics_pix_set_font_type(void *ctx, const char *type)
{
	context_t *c = ctx;
	c->font->type = type;
}

void nb_graphics_pix_set_font_size(void *ctx, uint16_t size)
{
	context_t *c = ctx;
	c->font->size = size;
}

void nb_graphics_pix_show_text(void *ctx, const char *str)
{
	/* PENDING */;
}

void nb_graphics_pix_get_text_attr(const void *ctx, const char *label,
				   nb_graphics_text_attr_t *attr)
{
	/* PENDING */;
}
