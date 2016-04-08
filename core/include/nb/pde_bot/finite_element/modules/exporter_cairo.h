#ifndef __NB_PDE_BOT_FINITE_ELEMENT_MODULES_EXPORTER_CAIRO_H__
#define __NB_PDE_BOT_FINITE_ELEMENT_MODULES_EXPORTER_CAIRO_H__

#include <stdint.h>
#include <stdbool.h>
#include "nb/geometric_bot.h"

void nb_fem_save_png(const vcn_msh3trg_t *const msh3trg,
		     const double *results,
		     const char* filename,
		     int width, int height);

#endif
