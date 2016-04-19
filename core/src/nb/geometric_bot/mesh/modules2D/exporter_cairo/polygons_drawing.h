#ifndef __NB_GEOMETRIC_BOT_MESH_MODULES2D_EXPORTER_CAIRO_POLYGONS_DRAW_H__
#define __NB_GEOMETRIC_BOT_MESH_MODULES2D_EXPORTER_CAIRO_POLYGONS_DRAW_H__

#include "nb/geometric_bot/mesh/elements2D/polygons.h"

void nb_mshpoly_export_png(const nb_mshpoly_t *const poly,
			   const char* filename, int width, int height);

void nb_mshpoly_export_eps(const nb_mshpoly_t *const poly,
			   const char* filename, int width, int height);

#endif
