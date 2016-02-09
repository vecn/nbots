/******************************************************************************
 *   Visual Cat: Visualization utilities for Cairo Graphics.                  *
 *   2011-2015 Victor Eduardo Cardoso Nungaray                                *
 *   Twitter: @victore_cardoso                                                *
 *   email: victorc@cimat.mx                                                  *
 ******************************************************************************/

/**
 * @file visual_cat-cairo.h
 * @brief Visualization utilities.
 * @author Victor Eduardo Cardoso Nungaray
 * @n victorc@@cimat.mx
 * @n <a href="https://twitter.com/victore_cardoso"> @@victore_cardoso </a>
 * @date November 2, 2015
 */

#ifndef __VCN_VISUAL_CAT_CAIRO_H__
#define __VCN_VISUAL_CAT_CAIRO_H__

#include <cairo.h>
#include "vcn/visual_cat.h"

void vcn_palette_draw_in_cairo(cairo_t *cr, 
			       const vcn_palette_t *const palette,
			       float x, float y, float w, float h,
			       float border, double min_v, double max_v);

#endif
