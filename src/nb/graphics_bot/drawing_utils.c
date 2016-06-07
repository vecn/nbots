#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "nb/graphics_bot/drawing_utils.h"

#include "nb/graphics_bot/drawing_tools.h"

#include "palette_struct.h"

#define MIN(a,b) (((a)<(b))?(a):(b))

static vcn_palette_t* palette_get_rainbow(void);
static vcn_palette_t* palette_get_sunset(void);
static vcn_palette_t* palette_get_french(void);

static void palette_draw_rectangle(void *draw_ptr,
				   const vcn_palette_t *const palette,
				   float x, float y, float w, float h,
				   float border);
static void palette_draw_zero_mark(void* draw_ptr,
				   const vcn_palette_t *const palette,
				   float x, float y, float w, float h,
				   double min_v, double max_v);
static void palette_draw_labels(void *draw_ptr, float font_size,
				float x, float y, float w, float h,
				double min_v, double max_v);

void nb_drawing_utils_set_center_and_zoom(camera_t *cam, const double box[4],
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

vcn_palette_t* vcn_palette_create(void)
{
	return calloc(1, sizeof(vcn_palette_t));
}

vcn_palette_t* vcn_palette_create_preset(nb_palette_preset preset)
{
	if (NB_RAINBOW == preset)
		return palette_get_rainbow();

	if (NB_SUNSET == preset)
		return palette_get_sunset();

	if (NB_FRENCH == preset)
		return palette_get_french();

	return NULL;
}

void vcn_palette_destroy(vcn_palette_t* palette)
{
	if (palette->ntics > 0) {
		free(palette->tics);
		free(palette->rgb);
	}
	free(palette);
}

void vcn_palette_clear(vcn_palette_t *palette)
{
	if (palette->ntics > 0) {
		free(palette->tics);
		free(palette->rgb);
		palette->tics = NULL;
		palette->rgb = NULL;
	}
	palette->ntics = 0;
}

void vcn_palette_add_colour(vcn_palette_t* palette, float tic,
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

void vcn_palette_get_colour(const vcn_palette_t *const palette,
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

static vcn_palette_t* palette_get_rainbow()
{
	vcn_palette_t* palette = vcn_palette_create();
	vcn_palette_add_colour(palette, 0.00f,   0,   0, 128);
	vcn_palette_add_colour(palette, 0.10f,   0,   0, 255);
	vcn_palette_add_colour(palette, 0.20f,   0, 128, 255);
	vcn_palette_add_colour(palette, 0.37f,   0, 255, 255);
	vcn_palette_add_colour(palette, 0.50f,   0, 255,   0);
	vcn_palette_add_colour(palette, 0.63f, 255, 255,   0);
	vcn_palette_add_colour(palette, 0.80f, 255, 128,   0);
	vcn_palette_add_colour(palette, 0.90f, 255,   0,   0);
	vcn_palette_add_colour(palette, 1.00f, 100,   0,   0);
	return palette;
}

static vcn_palette_t* palette_get_sunset()
{
	vcn_palette_t* palette = vcn_palette_create();
	vcn_palette_add_colour(palette, 0.00f,   0,   0,   0);
	vcn_palette_add_colour(palette, 0.15f,  20,   0, 100);
	vcn_palette_add_colour(palette, 0.30f, 100,   0, 200);
	vcn_palette_add_colour(palette, 0.80f, 220, 100,   0);
	vcn_palette_add_colour(palette, 1.00f, 255, 255,   0);
	return palette;
}

static vcn_palette_t* palette_get_french()
{
	vcn_palette_t* palette = vcn_palette_create();
	vcn_palette_add_colour(palette, 0.00f,   0,   0, 150);
	vcn_palette_add_colour(palette, 0.20f,   0,   0, 255);
	vcn_palette_add_colour(palette, 0.30f, 180, 180, 255);
	vcn_palette_add_colour(palette, 0.50f, 255, 255, 255);
	vcn_palette_add_colour(palette, 0.70f, 255, 180, 180);
	vcn_palette_add_colour(palette, 0.80f, 255,   0,   0);
	vcn_palette_add_colour(palette, 1.00f, 150,   0,   0);
	return palette;
}

void nb_drawing_draw_palette(void *draw_ptr, 
			     const vcn_palette_t *const palette,
			     float x, float y, float w, float h,
			     float border, double min_v, double max_v)
{
	palette_draw_rectangle(draw_ptr, palette, x, y, w, h, border);
	palette_draw_zero_mark(draw_ptr, palette, x, y, w, h, min_v, max_v);
	palette_draw_labels(draw_ptr, 10.0f, x, y, w, h, min_v, max_v);
}

static void palette_draw_rectangle(void *draw_ptr,
				   const vcn_palette_t *const palette,
				   float x, float y, float w, float h,
				   float border)
{
	nb_pattern_t *pat;
	pat = nb_pattern_create_linear(x, y+h, x, y);
	for (int i = 0; i < palette->ntics; i++) {
		nb_pattern_add_color_stop_rgb(pat, palette->tics[i],
					      palette->rgb[i * 3]/255.0,
					      palette->rgb[i*3+1]/255.0,
					      palette->rgb[i*3+2]/255.0);
	}
  
	if (0.0f < border) {
		nb_drawing_set_source_rgb(draw_ptr, 0.0, 0.0, 0.0);
		nb_drawing_raw_set_rectangle(draw_ptr, x-border, y-border,
					     w+2*border, h+2*border);
		nb_drawing_fill(draw_ptr);
	}
	nb_drawing_raw_set_rectangle(draw_ptr, x, y, w, h);
	nb_drawing_set_source(draw_ptr, pat);
	nb_drawing_fill(draw_ptr);
	nb_pattern_destroy(pat);
}

static void palette_draw_zero_mark(void* draw_ptr,
				   const vcn_palette_t *const palette,
				   float x, float y, float w, float h,
				   double min_v, double max_v)
{
	if (min_v * max_v >= 0.0)
		return;

	float factor = - min_v / (max_v - min_v);
	uint8_t rgb[3];
	vcn_palette_get_colour(palette, factor, rgb);
    
	float rcol = rgb[0] / 255.0f + 0.5f;
	float gcol = rgb[1] / 255.0f + 0.5f;
	float bcol = rgb[2] / 255.0f + 0.5f;

	if (rcol > 1.0)
		rcol -= 1.0;
	if (gcol > 1.0)
		gcol -= 1.0;
	if (bcol > 1.0)
		bcol -= 1.0;

	nb_drawing_set_source_rgb(draw_ptr, rcol, gcol, bcol);

	double yzero = h * factor;
	nb_drawing_raw_move_to(draw_ptr, x, y + h - yzero);
	nb_drawing_raw_line_to(draw_ptr, x + w, y + h - yzero);
	nb_drawing_stroke(draw_ptr);
}

static void palette_draw_labels(void *draw_ptr, float font_size,
				float x, float y, float w, float h,
				double min_v, double max_v)
{
	nb_drawing_set_font_type(draw_ptr, "Sans");
	nb_drawing_set_font_size(draw_ptr, font_size);
	nb_text_attr_t text_attr;;
	nb_drawing_set_source_rgb(draw_ptr, 0.0, 0.0, 0.0);

	char label[15];
	int n_labels = MIN((int)(h / (font_size + 2)), 10);
	double step_v = (max_v - min_v) / (n_labels - 1.0);
	for (int i = 0; i < n_labels; i++) {
		sprintf(label, "%.3e", max_v - i * step_v);
		nb_drawing_get_text_attr(draw_ptr, label, &text_attr);

		nb_drawing_raw_move_to(draw_ptr, x + w + 5.0f, 
				       y + text_attr.height/2.0 + 
				       i * h / (n_labels-1.0));
		nb_drawing_show_text(draw_ptr, label);
	}
}
