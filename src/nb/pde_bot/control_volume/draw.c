#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "nb/memory_bot.h"
#include "nb/graphics_bot.h"
#include "nb/geometric_bot.h"
#include "nb/pde_bot/control_volume/draw.h"

#include "integration_mesh.h"

static void draw(nb_graphics_context_t *g, int width, int height,
		 const void *draw_data);
static void set_camera(nb_graphics_context_t *g, int width, int height,
		       const nb_partition_t *part);

void nb_cvfa_draw_integration_mesh(const nb_partition_t *const part,
				   const char *filename, int w, int h)
{
	nb_graphics_export(filename, w, h, draw, part);
}

static void draw(nb_graphics_context_t *g, int width, int height,
		 const void *data_ptr)
{
	const nb_partition_t *part = data_ptr;

	uint32_t memsize = nb_cvfa_get_integration_mesh_memsize();
	char *memblock = nb_soft_allocate_mem(memsize);
	nb_partition_t *intmsh = (void*) memblock;

	nb_cvfa_init_integration_mesh(intmsh);
	nb_cvfa_load_integration_mesh(part, intmsh);

	if (!nb_graphics_is_camera_enabled(g))
		set_camera(g, width, height, part);

	nb_graphics_set_source(g, NB_LIGHT_BLUE);
	nb_partition_fill_elems(part, g);

	nb_graphics_set_source(g, NB_DARK_GRAY);
	nb_graphics_set_line_width(g, 1.0);
	nb_partition_draw_wires(part, g);

	nb_graphics_set_source(g, NB_BLUE);
	nb_graphics_set_line_width(g, 0.5);
	nb_partition_draw_wires(intmsh, g);

	nb_soft_free_mem(memsize, memblock);
}

static void set_camera(nb_graphics_context_t *g, int width, int height,
		       const nb_partition_t *part)
{
	double box[4];
	nb_partition_get_enveloping_box(part, box);

	nb_graphics_enable_camera(g);
	nb_graphics_camera_t* cam = nb_graphics_get_camera(g);
	nb_graphics_cam_fit_box(cam, box, width, height);
}
