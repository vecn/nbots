#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "nb/math_bot.h"
#include "nb/memory_bot.h"
#include "nb/graphics_bot/drawing_tools.h"

#include "eps_drawing.h"

typedef struct {
	FILE *fp; /* TEMPORAL FILE stream */	
} context_t;

static void set_eps_header(context_t *ctx);

void* nb_graphics_eps_create_context(int width, int height)
{
	uint32_t memsize = sizeof(context_t);
	char * memblock = nb_allocate_mem(memsize);
	context_t *ctx = (void*) memblock;
	ctx->fp = fopen("../../../TEMPORAL.eps", "w");
	set_eps_header(ctx);
	return ctx;
}

static void set_eps_header(context_t *c)
{
	fprintf(c->fp, "%%!PS-Adobe EPSF-3.0\n");
	fprintf(c->fp, "%%%%Creator: NBots -> Graphic Bot\n");
	fprintf(c->fp, "%%%%Title: NBots -> Graphic Bot\n");      /* TEMPORAL */
	fprintf(c->fp, "%%%%CreationDate: 04/06/2016 16:48:04\n");/* TEMPORAL */
	fprintf(c->fp, "%%%%DocumentData: Clean7Bit\n");          /* TEMPORAL */
	fprintf(c->fp, "%%%%Origin: 0 0\n");                      /* TEMPORAL */
	fprintf(c->fp, "%%%%BoundingBox: 0 0 1000 800\n");        /* TEMPORAL */
	fprintf(c->fp, "%%%%EndCommens\n");
}

void nb_graphics_eps_destroy_context(void *ctx)
{
	nb_free_mem(ctx);
}

void nb_graphics_eps_export_context(const void *ctx, const char *filename)
{
	const context_t *c = ctx;
	fprintf(c->fp, "showpage\n");
	fprintf(c->fp, "%%%%EOF\n");
	fclose(c->fp);

	/* TEMPORAL (COPY TO EXIT FILE) ****************************/
	FILE *source = fopen("../../../TEMPORAL.eps", "r");	 /**/
	FILE *target = fopen(filename, "w");			 /**/
	char ch;						 /**/
	while (( ch = fgetc(source) ) != EOF )			 /**/
		fputc(ch, target);				 /**/
	fclose(source);						 /**/
	fclose(target);						 /**/
	/***********************************************************/
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
	context_t *c = ctx;/* TEMPORAL BAD APPROX */
	fprintf(c->fp, "%g %g %g %g %g %g curveto\n",
		xcontrol, ycontrol, xcontrol, ycontrol, x, y);
}

void nb_graphics_eps_qrcurve_to(void *ctx, float x, float y,
				float xcontrol, float ycontrol, float w)
{
	context_t *c = ctx;/* TEMPORAL BAD APPROX */
	fprintf(c->fp, "%g %g %g %g %g %g curveto\n",
		xcontrol, ycontrol, xcontrol, ycontrol, x, y);
}

void nb_graphics_eps_curve_to(void *ctx, float x, float y,
			      float x0_control, float y0_control,
			      float x1_control, float y1_control)
{
	context_t *c = ctx;
	fprintf(c->fp, "%g %g %g %g %g %g curveto\n",
		x0_control, y0_control, x1_control, y1_control, x, y);
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
	context_t *c = ctx;
	fprintf(c->fp, "%g %g %g setrgbcolor\n",
		r / 255.0, g / 255.0, b / 255.0);
}

void nb_graphics_eps_set_source_rgba(void *ctx, uint8_t r, uint8_t g,
				     uint8_t b, uint8_t a)
{
	context_t *c = ctx;
	fprintf(c->fp, "%g %g %g setrgbcolor\n",
		r / 255.0, g / 255.0, b / 255.0);
}

void nb_graphics_eps_set_source_grad(void *ctx,
				     nb_graphics_grad_t grad,
				     float x1, float y1,
				     float x2, float y2,
				     nb_palette_t *pat)
{
	context_t *c = ctx;
	fprintf(c->fp, "%g %g %g setrgbcolor\n", 1.0f, 1.0f, 1.0f);
}

void nb_graphics_eps_set_source_trg(void *ctx,
				    float x1, float y1,
				    float x2, float y2,
				    float x3, float y3,
				    const uint8_t rgba1[4],
				    const uint8_t rgba2[4],
				    const uint8_t rgba3[4])
{
	context_t *c = ctx;
	float r = (rgba1[0] / 255.0 + rgba2[0] / 255.0 + rgba3[0] / 255.0);
	float g = (rgba1[1] / 255.0 + rgba2[1] / 255.0 + rgba3[1] / 255.0);
	float b = (rgba1[2] / 255.0 + rgba2[2] / 255.0 + rgba3[2] / 255.0);
	fprintf(c->fp, "%g %g %g setrgbcolor\n", r, g, b);
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
	context_t *c = ctx;
	fprintf(c->fp, "/%s findfont\n");
}

void nb_graphics_eps_set_font_size(void *ctx, uint16_t size)
{
	context_t *c = ctx;
	fprintf(c->fp, "%i scalefont\n", size);
}

void nb_graphics_eps_show_text(void *ctx, int x, int y, const char *str)
{
	context_t *c = ctx;
	fprintf(c->fp, "setfont\n");
	fprintf(c->fp, "(%s) show\n");
}

void nb_graphics_eps_get_text_attr(const void *ctx, const char *label,
				   nb_graphics_text_attr_t *attr)
{
	attr->width = strlen(label) * 10;
	attr->height = 10;/* TEMPORAL */
}
