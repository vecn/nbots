/******************************************************************************
 *   Visual Cat: Visualization utilities.                                     *
 *   2011-2015 Victor Eduardo Cardoso Nungaray                                *
 *   Twitter: @victore_cardoso                                                *
 *   email: victorc@cimat.mx                                                  *
 ******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "nb/math_cat.h"
#include "nb/visual_cat-cairo.h"

struct vcn_palette_s{
  /* The palette defines a serie of RGB colors to
   * colorize values in [0,1]
   *
   *  c1    c2       c3         c4  <- RGB colors
   *   |_____|________|__________|
   *   0    0.25     0.57        1  <- Tics
   */
  uchar ntics;   /* Number of tics */
  float* tics;   /* Sorted tics in [0,1] */
  uchar* rgb;    /* RGB Colors definition */
};

static void palette_draw_rectangle(cairo_t *cr,
				   const vcn_palette_t *const palette,
				   float x, float y, float w, float h,
				   float border);
static void palette_draw_zero_mark(cairo_t* cr,
				   const vcn_palette_t *const palette,
				   float x, float y, float w, float h,
				   double min_v, double max_v);
static void palette_draw_labels(cairo_t *cr, float font_size,
				float x, float y, float w, float h,
				double min_v, double max_v);

/******************* Implementation of public functions *******************/
void vcn_palette_draw_in_cairo(cairo_t* cr,
			       const vcn_palette_t *const palette,
			       float x, float y, float w, float h,
			       float border, double min_v, double max_v){
  palette_draw_rectangle(cr, palette, x, y, w, h, border);
  palette_draw_zero_mark(cr, palette, x, y, w, h, min_v, max_v);
  palette_draw_labels(cr, 10.0f, x, y, w, h, min_v, max_v);
}

/******************* Implementation of private functions *******************/
static void palette_draw_rectangle(cairo_t *cr,
				   const vcn_palette_t *const palette,
				   float x, float y, float w, float h,
				   float border)
{
  cairo_pattern_t *pat;
  pat = cairo_pattern_create_linear(x, y+h, x, y);
  for (int i = 0; i < palette->ntics; i++) {
    cairo_pattern_add_color_stop_rgb(pat, palette->tics[i],
				     palette->rgb[i * 3]/255.0,
				     palette->rgb[i*3+1]/255.0,
				     palette->rgb[i*3+2]/255.0);
  }
  
  if (0.0f < border) {
    cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
    cairo_rectangle(cr, x-border, y-border, w+2*border, h+2*border);
    cairo_fill(cr);
  }
  cairo_rectangle(cr, x, y, w, h);
  cairo_set_source(cr, pat);
  cairo_fill(cr);
  cairo_pattern_destroy(pat);
}

static void palette_draw_zero_mark(cairo_t* cr,
				   const vcn_palette_t *const palette,
				   float x, float y, float w, float h,
				   double min_v, double max_v)
{
  if (min_v * max_v >= 0.0)
    return;

  double factor = - min_v / (max_v - min_v);
  uchar rgb[3];
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

  cairo_set_source_rgb(cr, rcol, gcol, bcol);

  double yzero = h * factor;
  cairo_move_to(cr, x, y + h - yzero);
  cairo_line_to(cr, x + w, y + h - yzero);
  cairo_stroke(cr);
}

static void palette_draw_labels(cairo_t *cr, float font_size,
				float x, float y, float w, float h,
				double min_v, double max_v)
{
  cairo_select_font_face(cr, "Sans",
			 CAIRO_FONT_SLANT_NORMAL,
			 CAIRO_FONT_WEIGHT_NORMAL);
  cairo_set_font_size(cr, font_size);
  cairo_text_extents_t extents;
  cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);

  char label[15];
  int n_labels = min((int)(h / (font_size + 2)), 10);
  double step_v = (max_v - min_v) / (n_labels - 1.0);
  for (int i = 0; i < n_labels; i++) {
    sprintf(label, "%.3e", max_v - i * step_v);
    cairo_text_extents(cr, label, &extents);

    cairo_move_to(cr, x + w + 5.0f, 
		  y + extents.height/2.0 + 
		  i * h / (n_labels-1.0));
    cairo_show_text(cr, label);
  }
}
