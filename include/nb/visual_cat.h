/******************************************************************************
 *   Visual Cat: Visualization utilities.                                     *
 *   2011-2015 Victor Eduardo Cardoso Nungaray                                *
 *   Twitter: @victore_cardoso                                                *
 *   email: victorc@cimat.mx                                                  *
 ******************************************************************************/

/**
 * @file visual_cat.h
 * @brief Visualization utilities.
 * @author Victor Eduardo Cardoso Nungaray
 * @n victorc@@cimat.mx
 * @n <a href="https://twitter.com/victore_cardoso"> @@victore_cardoso </a>
 * @date September 9, 2015
 */

#ifndef __NB_VISUAL_CAT_H__
#define __NB_VISUAL_CAT_H__

enum {
	NB_PALETTE_RAINBOW,
	NB_PALETTE_SUNSET,
	NB_PALETTE_FRENCH
};

#include <stdint.h>

typedef struct vcn_palette_s vcn_palette_t;

vcn_palette_t* vcn_palette_create();
vcn_palette_t* vcn_palette_create_preset(int palette_id);
void vcn_palette_destroy(vcn_palette_t* palette);
void vcn_palette_clear(vcn_palette_t* palette);
void vcn_palette_add_colour(vcn_palette_t* palette, float tic,
			    uint8_t r, uint8_t g, uint8_t b);
void vcn_palette_get_colour(const vcn_palette_t *const palette,
			    float factor,
			    uint8_t rgb[3]);


#endif
