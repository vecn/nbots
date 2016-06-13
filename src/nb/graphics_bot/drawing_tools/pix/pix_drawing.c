#include <stdlib.h>
#include <stdint.h>
#include <alloca.h>

#include "nb/math_bot.h"
#include "nb/container_bot.h"
#include "nb/image_bot.h"
#include "nb/graphics_bot/drawing_tools.h"

#include "pix_drawing.h"

#define TURTLE_STRUCT NB_QUEUE

typedef struct {
	vcn_image_t *img;
	nb_container_t *turtle;
} context_t;

void* nb_graphics_pix_create_context(int width, int height)
{
	uint32_t ctx_size = sizeof(context_t);
	uint32_t img_size = sizeof(vcn_image_t);
	uint32_t pix_size = 4 * width * height;
	uint32_t cnt_size = nb_container_get_memsize(TURTLE_STRUCT);
	uint32_t memsize = ctx_size + img_size + pix_size + cnt_size;
	char *memblock = malloc(memsize);
	context_t *ctx = (void*) memblock;

	ctx->img = (void*) (memblock + ctx_size);
	vcn_image_init(ctx->img);
	ctx->img->width = width;
	ctx->img->height = height;
	ctx->img->comp_x_pixel = 4;
	ctx->img->pixels = (void*) (memblock + ctx_size + img_size);

	ctx->turtle = (void*) (memblock + ctx_size + img_size + pix_size);
	nb_container_init(ctx->turtle, TURTLE_STRUCT);
	
	return ctx;
}

void nb_graphics_pix_destroy_context(void *ctx)
{
	free(ctx);
}

void nb_graphics_pix_export_context(void *ctx, const char *filename)
{
	/* PNG, BMP, TXT */;
}

void nb_graphics_pix_move_to(void *ctx, double x, double y)
{
	;
}

void nb_graphics_pix_line_to(void *ctx, double x, double y)
{
	;
}

void nb_graphics_pix_arc(void *ctx,
			 double x, double y, double r,
			 double a0, double a1)
{
	;
}

void nb_graphics_pix_set_circle(void *ctx,
				double x, double y, double r)
{

}

void nb_graphics_pix_set_rectangle(void *ctx, double x1, double y1,
				       double x2, double y2)
{
	;
}

void nb_graphics_pix_close_path(void *ctx)
{
	;
}

void nb_graphics_pix_set_line_width(void *ctx, double w)
{
	;
}

void nb_graphics_pix_set_source_rgb(void *ctx,
				    uint8_t r, uint8_t g, uint8_t b)
{
	;
}

void nb_graphics_pix_set_source_rgba(void *ctx, uint8_t r, uint8_t g,
				     uint8_t b, uint8_t a)
{
	;
}

void nb_graphics_pix_set_source_grad(void *ctx,
				     nb_graphics_grad_t grad,
				     double x1, double y1,
				     double x2, double y2,
				     nb_graphics_palette_t *pat)
{
	;
}

void nb_graphics_pix_set_source_trg(void *ctx,
				    double x1, double y1,
				    double x2, double y2,
				    double x3, double y3,
				    uint8_t r1, uint8_t g1, uint8_t b1,
				    uint8_t r2, uint8_t g2, uint8_t b2,
				    uint8_t r3, uint8_t g3, uint8_t b3)
{
	;
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

void nb_graphics_pix_get_text_attr(void *ctx, const char *label,
				   nb_text_attr_t *attr)
{
	;
}
