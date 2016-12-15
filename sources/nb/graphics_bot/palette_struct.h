#ifndef __NB_GRAPHICS_BOT_PALETTE_STRUCT_H__
#define __NB_GRAPHICS_BOT_PALETTE_STRUCT_H__

#include <stdint.h>

struct nb_palette_s {
	/* The palette defines a serie of RGB colors to
	 * colorize values in [0,1]
	 *
	 *  c1    c2       c3         c4  <- RGB colors
	 *   |_____|________|__________|
	 *   0    0.25     0.57        1  <- Tics
	 */
	uint8_t ntics;       /* Number of tics */
	float tics[10];      /* Sorted tics in [0,1] */
	uint8_t rgba[40];    /* RGBA Colors definition */
};

#endif
