#include <stdint.h>

#include "nb/math_bot.h"
#include "nb/memory_bot.h"
#include "nb/graphics_bot/drawing_tools.h"

#include "eps_drawing.h"

enum {
	FILE *fp; /* TEMPORAL FILE stream */	
} context_t;

static void set_eps_header(context_t *ctx);

void* nb_graphics_eps_create_context(int width, int height)
{
	uint32_t memsize = sizeof(context_t);
	char * memblock = nb_allocate_mem(memsize);
	context_t *ctx = (void*) memblock;
	ctx->fp = fopen("../../../TEMPORAL.eps");
	set_eps_header(ctx);
	return ctx;
}

static void set_eps_header(context_t *c)
{
	fprintf(c->fp, "%!PS-Adobe EPSF-3.0\n");
	fprintf(c->fp, "%%Creator: NBots -> Graphic Bot\n");
	fprintf(c->fp, "%%Title: NBots -> Graphic Bot\n");      /* TEMPORAL */
	fprintf(c->fp, "%%CreationDate: 04/06/2016 16:48:04\n");/* TEMPORAL */
	fprintf(c->fp, "%%DocumentData: Clean7Bit\n");/* TEMPORAL */
	fprintf(c->fp, "%%Origin: 0 0\n");/* TEMPORAL */
	fprintf(c->fp, "%%BoundingBox: 0 0 400 400\n");         /* TEMPORAL */
	fprintf(c->fp, "%%EndCommens\n");
}

void nb_graphics_eps_destroy_context(void *ctx)
{
	free(ctx);
}

void nb_graphics_eps_export_context(const void *ctx, const char *filename)
{
	const context_t *c = ctx;
	fprintf(c->fp, "%%EOF\n");
	fclose(c->fp);
}

void nb_graphics_eps_move_to(void *ctx, float x, float y)
{
	context_t *c = ctx;
	fprintf(c->fp, "newpath\n");
	fprintf(c->fp, "%g %g moveto\n", x, y);
}

void nb_graphics_eps_line_to(void *ctx, float x, float y)
{
	context_t *c = ctx;
	fprintf(c->fp, "%g %g lineto\n", x, y);
}

void nb_graphics_eps_qcurve_to(void *ctx, float x, float y,
			       float xcontrol, float ycontrol)
{
	;
}

void nb_graphics_eps_qrcurve_to(void *ctx, float x, float y,
				float xcontrol, float ycontrol, float w)
{
	;
}

void nb_graphics_eps_curve_to(void *ctx, float x, float y,
			      float x0_control, float y0_control,
			      float x1_control, float y1_control)
{
	;
}

void nb_graphics_eps_close_path(void *ctx)
{
	context_t *c = ctx;
	fprintf(c->fp, "closepath\n");
}

void nb_graphics_eps_set_line_width(void *ctx, float w)
{
	context_t *c = ctx;
	fprintf(c->fp, "%g setlinewidth\n", w);
}

void nb_graphics_eps_set_source_rgb(void *ctx,
				    uint8_t r, uint8_t g, uint8_t b)
{
	;
}

void nb_graphics_eps_set_source_rgba(void *ctx, uint8_t r, uint8_t g,
				     uint8_t b, uint8_t a)
{
	;
}

void nb_graphics_eps_set_source_grad(void *ctx,
				     nb_graphics_grad_t grad,
				     float x1, float y1,
				     float x2, float y2,
				     nb_palette_t *pat)
{
	;
}

void nb_graphics_eps_set_source_trg(void *ctx,
				    float x1, float y1,
				    float x2, float y2,
				    float x3, float y3,
				    const uint8_t rgba1[4],
				    const uint8_t rgba2[4],
				    const uint8_t rgba3[4])
{
	;
}

void nb_graphics_eps_fill(void *ctx)
{
	context_t *c = ctx;
	fprintf(c->fp, "fill\n");
}

void nb_graphics_eps_fill_preserve(void *ctx)
{
	context_t *c = ctx;
	fprintf(c->fp, "gsave\n");
	fprintf(c->fp, "fill\n");
	fprintf(c->fp, "grestore\n");
}

void nb_graphics_eps_stroke(void *ctx)
{
	context_t *c = ctx;
	fprintf(c->fp, "stroke\n");
}

void nb_graphics_eps_stroke_preserve(void *ctx)
{
	context_t *c = ctx;
	fprintf(c->fp, "gsave\n");
	fprintf(c->fp, "stroke\n");
	fprintf(c->fp, "grestore\n");
}

void nb_graphics_eps_set_font_type(void *ctx, const char *type)
{
	;
}

void nb_graphics_eps_set_font_size(void *ctx, uint16_t size)
{
	;
}

void nb_graphics_eps_show_text(void *ctx, int x, int y, const char *str)
{
	;
}

void nb_graphics_eps_get_text_attr(const void *ctx, const char *label,
				   nb_graphics_text_attr_t *attr)
{
	;
}
