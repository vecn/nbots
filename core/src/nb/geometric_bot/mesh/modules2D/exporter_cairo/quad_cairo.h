#ifndef __NB_GEOMETRIC_BOT_MESH_MODULES2D_EXPORTER_CAIRO_QUAD_CAIRO_H__
#define __NB_GEOMETRIC_BOT_MESH_MODULES2D_EXPORTER_CAIRO_QUAD_CAIRO_H__

#include "nb/geometric_bot/mesh/elements2D/quad.h"

void nb_mshquad_export_png(const nb_mshquad_t *const quad,
			   const char* filename, int width, int height);

#endif
