#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <alloca.h>

#include "nb/math_bot.h"
#include "nb/image_bot.h"
#include "nb/graphics_bot/drawing_tools.h"

#include "pix_drawing.h"

#define TURTLE_STATIC_MEMSIZE 10
#define TURTLE_DYNAMIC_MEMINCREASE 15

enum {
	RGB, RGBA, GRAD, TRG
};

enum {
	MOVE_TO, LINE_TO, ARC, CLOSE_PATH
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
	nb_graphics_grad_t type;
	nb_graphics_palette_t *pal;
} source_grad_t;

typedef union {
	uint8_t comp[12];
	source_grad_t grad;
} source_color_t;

typedef struct {
	uint8_t type;
	float vtx[6];
	source_color_t color;
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
			   uint8_t r1, uint8_t g1, uint8_t b1,
			   uint8_t r2, uint8_t g2, uint8_t b2,
			   uint8_t r3, uint8_t g3, uint8_t b3);

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
	source->type = RGB;
	source->color.comp[0] = r;
	source->color.comp[1] = g;
	source->color.comp[2] = b;
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
	source_set_rgb(ctx->source, r, g, b);
}

void nb_graphics_pix_set_source_rgba(void *ctx, uint8_t r, uint8_t g,
				     uint8_t b, uint8_t a)
{
	context_t *c = ctx;
	source_set_rgba(ctx->source, r, g, b, a);
}

static void source_set_rgba(source_t *source, uint8_t r,
			    uint8_t g, uint8_t b, uint8_t a)
{
	memset(source, 0, sizeof(source_t));
	source->type = RGBA;
	source->color.comp[0] = r;
	source->color.comp[1] = g;
	source->color.comp[2] = b;
	source->color.comp[3] = a;
}
void nb_graphics_pix_set_source_grad(void *ctx,
				     nb_graphics_grad_t grad,
				     float x1, float y1,
				     float x2, float y2,
				     nb_graphics_palette_t *pal)
{
	context_t *c = ctx;
	source_set_grad(ctx->source, x1, y1, x2, y2, pal);
}

static void source_set_grad(source_t *source,
			    nb_graphics_grad_t grad,
			    float x1, float y1,
			    float x2, float y2,
			    nb_graphics_palette_t *pal)
{
	memset(source, 0, sizeof(source_t));
	source->type = GRAD;
	source->vtx[0] = x1;
	source->vtx[1] = y1;
	source->vtx[2] = x2;
	source->vtx[3] = y2;
	source->color.grad.pal = pal;
	source->color.grad.type = grad;
}

void nb_graphics_pix_set_source_trg(void *ctx,
				    float x1, float y1,
				    float x2, float y2,
				    float x3, float y3,
				    uint8_t r1, uint8_t g1, uint8_t b1,
				    uint8_t r2, uint8_t g2, uint8_t b2,
				    uint8_t r3, uint8_t g3, uint8_t b3)
{
	context_t *c = ctx;
	source_set_trg(ctx->source, x1, y1, x2, y2, x3, y3,
		       r1, g1, b1, r2, g2, b2, r3, g3, b3);
}

static void source_set_trg(source_t *source,
			   float x1, float y1,
			   float x2, float y2,
			   float x3, float y3,
			   uint8_t r1, uint8_t g1, uint8_t b1,
			   uint8_t r2, uint8_t g2, uint8_t b2,
			   uint8_t r3, uint8_t g3, uint8_t b3)
{
	memset(source, 0, sizeof(source_t));
	source->type = TRG;
	source->vtx[0] = x1;
	source->vtx[1] = y1;
	source->vtx[2] = x2;
	source->vtx[3] = y2;
	source->vtx[4] = x3;
	source->vtx[5] = y3;
	source->color.comp[0] = r1;
	source->color.comp[1] = g1;
	source->color.comp[2] = b1;
	source->color.comp[3] = 1;
	source->color.comp[4] = r2;
	source->color.comp[5] = g2;
	source->color.comp[6] = b2;
	source->color.comp[7] = 1;
	source->color.comp[8] = r3;
	source->color.comp[9] = g3;
	source->color.comp[10] = b3;
	source->color.comp[11] = 1;
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
	;
}

void nb_graphics_pix_set_font_type(void *ctx, const char *type)
{
	;
}

void nb_graphics_pix_set_font_size(void *ctx, uint16_t size)
{
	;
}

void nb_graphics_pix_show_text(void *ctx, const char *str)
{
	;
}

void nb_graphics_pix_get_text_attr(const void *ctx, const char *label,
				   nb_text_attr_t *attr)
{
	;
}
