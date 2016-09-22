#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>

#include "nb/geometric_bot.h"
#include "nb/graphics_bot.h"

static void draw_disks(const nb_mshpack_t *mshpack,
		       nb_graphics_context_t *g,
		       bool fill_disks,
		       void *source_data,
		       void (*set_source)(const nb_mshpack_t *msh,
					  nb_graphics_context_t *g,
					  uint32_t i, void *data));

static void set_source_null(const nb_mshpack_t *msh,
			    nb_graphics_context_t *g,
			    uint32_t i, void *data);

static void set_source_field(const nb_msh3trg_t *msh,
			     nb_graphics_context_t *g,
			     uint32_t i, void *data);

static void set_source_classes(const nb_msh3trg_t *msh,
			       nb_graphics_context_t *g,
			       uint32_t i, void *data);

void nb_mshpack_draw_wires(const void *msh,
			   nb_graphics_context_t *g)
{
	draw_disks(msh, g, false, NULL, NULL);
}

static void draw_disks(const nb_mshpack_t *mshpack,
		       nb_graphics_context_t *g,
		       bool fill_disks,
		       void *source_data,
		       void (*set_source)(const nb_mshpack_t *msh,
					  nb_graphics_context_t *g,
					  uint32_t i, void *data))
{
	uint32_t N_elems = nb_mshpack_get_N_elems(mshpack);
	for (uint32_t i = 0; i < N_elems; i++) {
		double x = nb_mshpack_elem_get_x(mshpack, i);
		double y = nb_mshpack_elem_get_y(mshpack, i);
		double r = nb_mshpack_elem_get_radius(mshpack, i);
		nb_graphics_set_circle(g, x, y, r);
		
		if (fill_disks) {
			set_source(mshpack, g, i, source_data);
			nb_graphics_fill(g);
		} else {
			nb_graphics_stroke(g);
		}
	}
}

void nb_mshpack_draw_boundaries(const void *msh,
				nb_graphics_context_t *g)
{
	;/* NULL statement */
}

void nb_mshpack_fill_elems(const void *msh,
			   nb_graphics_context_t *g)
{
	draw_disks(msh, g, true, NULL, set_source_null);
}

static void set_source_null(const nb_mshpack_t *msh,
			    nb_graphics_context_t *g,
			    uint32_t i, void *data)
{
	;/* NULL statement */
}

void nb_mshpack_fill_elems_field_on_nodes(const void *msh,
					  nb_graphics_context_t *g,
					  const double *normalized_field,
					  nb_graphics_palette_preset palette)
{
	;/* NULL statement */
}

void nb_mshpack_fill_elems_field_on_elems(const void *msh,
					  nb_graphics_context_t *g,
					  const double *normalized_field,
					  nb_graphics_palette_preset palette)
{
	const void *data[2];
	data[0] = (void*) normalized_field;
	data[1] = nb_graphics_palette_create_preset(palette);
	
	draw_disks(msh, g, true, NULL, set_source_field);
	
	nb_graphics_palette_destroy(data[1]);
}

static void set_source_field(const nb_msh3trg_t *msh,
			     nb_graphics_context_t *g,
			     uint32_t i, void *data)
{
	void **cls_data = data;
	double *field = cls_data[0];
	nb_graphics_palette_t *palette = cls_data[1];
	
	uint8_t c[4];
	nb_graphics_palette_get_rgba(palette, field[i], c);
	nb_graphics_set_source_rgba(g, c[0], c[1], c[2], c[3]);
}

void nb_mshpack_fill_elems_classes(const void *msh,
				   nb_graphics_context_t *g,
				   const uint8_t *class, uint8_t N_colors,
				   const nb_graphics_color_t *colors)
{
	void *data[3];
	data[0] = (void*)class;
	data[1] = (void*)colors;
	data[2] = &N_colors;
	draw_disks(msh, g, true, data, set_source_classes);
}

static void set_source_classes(const nb_msh3trg_t *msh,
			       nb_graphics_context_t *g,
			       uint32_t i, void *data)
{
	void **cls_data = data;
	uint8_t *class = cls_data[0];
	nb_graphics_color_t *colors = cls_data[1];
	uint8_t *N_colors = cls_data[2];
	
	uint8_t id_class = class[i];
	nb_graphics_color_t c = colors[id_class % *N_colors];

	nb_graphics_set_source(g, c);
}


void nb_mshpack_fill_nodes(const void *msh,
			   nb_graphics_context_t *g)
{
	;/* NULL statement */
}

void nb_mshpack_fill_nodes_classes(const void *msh,
				   nb_graphics_context_t *g,
				   const uint8_t *class, uint8_t N_colors,
				   const nb_graphics_color_t *colors)
{
	;/* NULL statement */
}
