#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <alloca.h>

#include "nb/math_bot.h"
#include "nb/image_bot.h"
#include "nb/graphics_bot/drawing_tools.h"

#include "../../palette_struct.h"

#include "rasterizer.h"

#include "pix_drawing.h"

#define TURTLE_STATIC_MEMSIZE 10       /* turtle_step */
#define TURTLE_DYNAMIC_MEMINCREASE 15  /* turtle_step */
#define PIXMASK_STACK_MEMSIZE 250      /* Bytes */

enum {
	SOLID, GRAD, TRG
};

enum {
	MOVE_TO, LINE_TO, ARC, CLOSE_PATH
	/* QCURVE_TO QRCURVE_TO CURVE_TO, ellipse() [TEMPORAL] */
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
} turtle_struct;

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
	vcn_image_t *img;
	turtle_struct *turtle;
	float line_width;
	source_t *source;
	font_t *font;
} context_t;

static void turtle_clear(turtle_struct *turtle);
static void turtle_add(turtle_struct *turtle, uint8_t type,
		       float x, float y, float v1, float v2,
		       float v3, float v4);
static void source_set_rgb(source_t *source, uint8_t r,
			   uint8_t g, uint8_t b);
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

void* nb_graphics_pix_create_context(int width, int height)
{
	uint32_t ctx_size = sizeof(context_t);
	uint32_t img_size = sizeof(vcn_image_t);
	uint32_t pix_size = 4 * width * height;
	uint32_t trt_size = sizeof(turtle_struct);
	uint16_t src_size = sizeof(source_t);
	uint16_t fnt_size = sizeof(font_t);
	uint32_t memsize = ctx_size + img_size + pix_size +
		trt_size + src_size + fnt_size;
	char *memblock = malloc(memsize);
	context_t *ctx = (void*) memblock;

	ctx->img = (void*) (memblock + ctx_size);
	vcn_image_init(ctx->img);
	ctx->img->width = width;
	ctx->img->height = height;
	ctx->img->comp_x_pixel = 4;
	ctx->img->pixels = (void*) (memblock + ctx_size + img_size);

	ctx->turtle = (void*) (memblock + ctx_size + img_size + pix_size);
	memset(ctx->turtle, 0, trt_size);

	ctx->line_width = 1.0;
	
	ctx->source = (void*) (memblock + ctx_size + img_size +
			       pix_size + trt_size);
	source_set_rgb(ctx->source, 0, 0, 0);
	
	ctx->font = (void*) (memblock + ctx_size + img_size +
			       pix_size + trt_size + src_size);
	ctx->font->type = "Mono";
	ctx->font->size = 18;

	return ctx;
}

static void source_set_rgb(source_t *source, uint8_t r,
			   uint8_t g, uint8_t b)
{
	memset(source, 0, sizeof(source_t));
	source->source_type = SOLID;
	source->pal = &(source->color);
	nb_graphics_palette_add_rgba(source->pal, 0.0,
				     r, g, b, 255);
}

void nb_graphics_pix_destroy_context(void *ctx)
{
	context_t *c = ctx;
	vcn_image_finish(c->img);
	turtle_clear(c->turtle);
	free(ctx);
}

static void turtle_clear(turtle_struct *turtle)
{
	if (turtle->N_alloc > 0)
		free(turtle->dynamic_mem);
	memset(turtle, 0, sizeof(turtle_struct));
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

static void turtle_add(turtle_struct *turtle, uint8_t type,
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

void nb_graphics_pix_arc(void *ctx,
			 float x, float y, float r,
			 float a0, float a1)
{
	context_t *c = ctx;
	turtle_add(c->turtle, ARC, x, y, r, a0, a1, 0);
}

void nb_graphics_pix_close_path(void *ctx)
{
	context_t *c = ctx;
	turtle_add(c->turtle, CLOSE_PATH,
		   0, 0, 0, 0, 0, 0);
}

void nb_graphics_pix_set_circle(void *ctx,
				float x, float y, float r)
{
	context_t *c = ctx;
	turtle_add(c->turtle, ARC, x, y, r,
		   0.0, 2 * NB_PI, 0);
}

void nb_graphics_pix_set_rectangle(void *ctx, float x1, float y1,
				   float x2, float y2)
{
	context_t *c = ctx;
	turtle_add(c->turtle, MOVE_TO, x1, y1,
		   0, 0, 0, 0);
	turtle_add(c->turtle, LINE_TO, x2, y1,
		   0, 0, 0, 0);
	turtle_add(c->turtle, LINE_TO, x2, y2,
		   0, 0, 0, 0);
	turtle_add(c->turtle, LINE_TO, x1, y2,
		   0, 0, 0, 0);
	turtle_add(c->turtle, CLOSE_PATH,
		   0, 0, 0, 0, 0, 0);
}

void nb_graphics_pix_set_line_width(void *ctx, float w)
{
	context_t *c = ctx;
	c->line_width = w;
}

void nb_graphics_pix_set_source_rgb(void *ctx,
				    uint8_t r, uint8_t g, uint8_t b)
{
	context_t *c = ctx;
	source_set_rgb(c->source, r, g, b);
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
	source->vtx[2] = x2;
	source->vtx[3] = y2;
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
	;
}

void nb_graphics_pix_fill_preserve(void *ctx)
{
	;
}

void nb_graphics_pix_stroke(void *ctx)
{
	context_t *c = ctx;
	/* PENDING */
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
	;
}

void nb_graphics_pix_get_text_attr(const void *ctx, const char *label,
				   nb_graphics_text_attr_t *attr)
{
	;
}
