#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>

#include "nb/math_bot.h"
#include "nb/container_bot.h"
#include "nb/geometric_bot.h"

#include "drawing_utils.h"
#include "drawing_tools.h"

#include "../../mesh2D_structs.h"

#include "polygons_drawing.h"

static void draw_mesh(void *draw_ptr, int width, int height,
		      const void *const poly_ptr);

void nb_mshpoly_export_png(const nb_mshpoly_t *const poly,
			   const char* filename, int width, int height)
{
	nb_drawing_export_png(filename, width, height, draw_mesh, poly);
}

void nb_mshpoly_export_eps(const nb_mshpoly_t *const poly,
			   const char* filename, int width, int height)
{
	nb_drawing_export_eps(filename, width, height, draw_mesh, poly);
}


static void draw_mesh(void *draw_ptr, int width, int height,
		      const void *const poly_ptr)
{
	;
}
