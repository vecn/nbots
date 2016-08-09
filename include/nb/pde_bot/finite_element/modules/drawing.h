#ifndef __NB_PDE_BOT_FINITE_ELEMENT_MODULES_DRAWING_H__
#define __NB_PDE_BOT_FINITE_ELEMENT_MODULES_DRAWING_H__

#include <stdint.h>
#include <stdbool.h>
#include "nb/geometric_bot.h"

void nb_fem_save(const vcn_msh3trg_t *const msh3trg,
		 const double *distortion, /* Can be NULL */
		 double max_distortion,
		 const double *results,
		 const char* filename,
		 int width, int height);

#endif
