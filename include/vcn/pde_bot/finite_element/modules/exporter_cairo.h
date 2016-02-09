#ifndef __NBOTS_PDE_BOT_FINITE_ELEMENT_MODULES_EXPORTER_CAIRO_H__
#define __NBOTS_PDE_BOT_FINITE_ELEMENT_MODULES_EXPORTER_CAIRO_H__

#include <stdint.h>
#include <stdbool.h>
#include "vcn/geometric_bot.h"

void nb_fem_save_png(const vcn_mesh_t *const mesh,
		     const double *results,
		     const char* filename,
		     int width, int height);

void nb_fem_save_eps(const vcn_mesh_t *const mesh,
		     const double *results,
		     const char* filename,
		     int width, int height);

#endif
