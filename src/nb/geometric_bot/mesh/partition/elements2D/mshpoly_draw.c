#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>

#include "nb/math_bot.h"
#include "nb/container_bot.h"
#include "nb/geometric_bot.h"
#include "nb/graphics_bot.h"

static void draw_sgm(const nb_mshpoly_t *poly,
		     nb_graphics_context_t *g,
		     uint32_t isgm);
static void set_source_null(const nb_mshpoly_t *msh,
			    nb_graphics_context_t *g,
			    uint32_t i, void *data);
static void fill_elems(const nb_mshpoly_t *poly,
		       nb_graphics_context_t *g,
		       const void *source_data,
		       void (*set_source)(const nb_mshpoly_t *msh,
					  nb_graphics_context_t *g,
					  uint32_t i, void *data));
static void set_source_field_on_nodes(const nb_mshpoly_t *msh,
				      nb_graphics_context_t *g,
				      uint32_t i, void *data);
static void set_source_field_on_elems(const nb_mshpoly_t *msh,
				      nb_graphics_context_t *g,
				      uint32_t i, void *data);
static void set_source_classes(const nb_mshpoly_t *msh,
			       nb_graphics_context_t *g,
			       uint32_t i, void *data);
static void fill_nodes(const nb_mshpoly_t *msh,
		       nb_graphics_context_t *g,
		       const void *source_data,
		       void (*set_source)(const nb_mshpoly_t *msh,
					  nb_graphics_context_t *g,
					  uint32_t i, void *data));

void nb_mshpoly_draw_wires(const void *msh,
			   nb_graphics_context_t *g)
{
	uint32_t N_edg = nb_msh3trg_get_N_edges(msh);
	for (uint32_t i = 0; i < N_edg; i++) {
		uint32_t n1 = nb_msh3trg_edge_get_1n(msh, i);
		uint32_t n2 = nb_msh3trg_edge_get_2n(msh, i);

		double x = nb_msh3trg_node_get_x(msh, n1);
		double y = nb_msh3trg_node_get_y(msh, n1);
		nb_graphics_move_to(g, x, y);

		x = nb_msh3trg_node_get_x(msh, n2);
		y = nb_msh3trg_node_get_y(msh, n2);
		nb_graphics_line_to(g, x, y);

		nb_graphics_stroke(g);
	}
}

void nb_mshpoly_draw_boundaries(const void *msh,
				nb_graphics_context_t *g)
{
	uint32_t N_sgm = nb_mshpoly_get_N_insgm(msh);
	for (uint32_t i = 0; i < N_sgm; i++)
		draw_sgm(msh, g, i);
}

static void draw_sgm(const nb_mshpoly_t *poly,
		     nb_graphics_context_t *g,
		     uint32_t isgm)
{
	uint32_t N_vtx = nb_mshpoly_insgm_get_N_nodes(poly, isgm);
	if (0 < N_vtx) {
		uint32_t nid = nb_mshpoly_get_node_x_insgm(poly, isgm, 0);
		double x = nb_mshpoly_get_x_node(poly, nid);
		double y = nb_mshpoly_get_y_node(poly, nid);
		nb_graphics_move_to(g, x, y);

		for (uint32_t j = 1; j < N_vtx; j++) {
			nid = nb_mshpoly_get_node_x_insgm(poly, isgm, j);
			x = nb_mshpoly_get_x_node(poly, nj);
			y = nb_mshpoly_get_y_node(poly, nj);
			nb_graphics_line_to(g, x, y);
		}
		nb_graphics_stroke(g);
	}
}

void nb_mshpoly_fill_elems(const void *msh,
			   nb_graphics_context_t *g)
{
	fill_elems(msh, g, NULL, set_source_null);
}

static void set_source_null(const nb_mshpoly_t *msh,
			    nb_graphics_context_t *g,
			    uint32_t i, void *data)
{
	;/* NULL statement */
}


static void fill_elems(const nb_mshpoly_t *poly,
		       nb_graphics_context_t *g,
		       const void *source_data,
		       void (*set_source)(const nb_mshpoly_t *msh,
					  nb_graphics_context_t *g,
					  uint32_t i, void *data))
{
	uint32_t N_elems = nb_mshpoly_get_N_elems(poly);
	for (uint32_t i = 0; i < N_elems; i++) {
		uint32_t id = nb_mshpoly_elem_get_adj(poly, i, 0);
		double x = nb_mshpoly_get_x_node(poly, id);
		double y = nb_mshpoly_get_y_node(poly, id);
		nb_graphics_move_to(g, x, y);

		uint16_t N_adj = nb_mshpoly_elem_get_N_adj(poly, i);
		for (uint16_t j = 1; j < N_adj; j++) {
			id = nb_mshpoly_elem_get_adj(poly, i, j);
			x = nb_mshpoly_get_x_node(poly, id);
			y = nb_mshpoly_get_y_node(poly, id);
			nb_graphics_line_to(g, x, y);			
		}
		nb_graphics_close_path(g);

		set_source(poly, g, i, source_data);

		nb_graphics_fill(g);
	}
}

void nb_mshpoly_fill_elems_field_on_nodes(const void *msh,
					  nb_graphics_context_t *g,
					  const double *normalized_field,
					  nb_graphics_palette_preset palette)
{
	const void *data[2];
	data[0] = (void*) normalized_field;
	data[1] = nb_graphics_palette_create_preset(palette);
	
	fill_elems(msh, g, data, set_source_field_on_nodes);
	
	nb_graphics_palette_destroy(data[1]);
}

static void set_source_field_on_nodes(const nb_mshpoly_t *msh,
				      nb_graphics_context_t *g,
				      uint32_t i, void *data)
{
	void **cls_data = data;
	double *field = data[0];
	nb_graphics_palette_t *palette = data[1];

	uint16_t N_adj = nb_mshpoly_elem_get_N_adj(msh, i);
	double *vtx = alloca(2 * N_adj);
	uint8_t *rgba = alloca(4 * N_adj);

	for (uint32_t i = 0; i < N_adj; i++) {
		uint32_t n = nb_msh3trg_elem_get_adj(msh, i, j);
		vtx[i * 2] = nb_msh3trg_node_get_x(msh, n);
		vtx[i*2+1] = nb_msh3trg_node_get_y(msh, n);
		nb_graphics_palette_get_rgba(palette, field[n],
					     &(rgba[i*4]));
	}

	nb_graphics_set_source_fan(g, N, vtx, rgba);
}

void nb_mshpoly_fill_elems_field_on_elems(const void *msh,
					  nb_graphics_context_t *g,
					  const double *normalized_field,
					  nb_graphics_palette_preset palette)
{
	const void *data[2];
	data[0] = (void*) normalized_field;
	data[1] = nb_graphics_palette_create_preset(palette);
	
	fill_elems(msh, g, data, set_source_field_on_elems);
	
	nb_graphics_palette_destroy(data[1]);
}

static void set_source_field_on_elems(const nb_mshpoly_t *msh,
				      nb_graphics_context_t *g,
				      uint32_t i, void *data)
{
	void **cls_data = data;
	double *field = data[0];
	nb_graphics_palette_t *palette = data[1];
	
	uint8_t c[4];
	nb_graphics_palette_get_rgba(palette, field[i], c);
	nb_graphics_set_source_rgba(g, c[0], c[1], c[2], c[3]);
}

void nb_mshpoly_fill_elems_classes(const void *msh,
				   nb_graphics_context_t *g,
				   const uint8_t *class, uint8_t N_colors,
				   const nb_graphics_color_t *colors)
{
	void *data[3];
	data[0] = (void*)class;
	data[1] = (void*)colors;
	data[2] = &N_colors;
	fill_elems(msh, g, data, set_source_classes);
}

static void set_source_classes(const nb_mshpoly_t *msh,
			       nb_graphics_context_t *g,
			       uint32_t i, void *data)
{
	void **cls_data = data;
	uint8_t *class = cls_data[0];
	nb_graphics_color_t *colors = cls_data[1];
	uint8_t *N_colors = cls_data[2];
	
	uint8_t id_class = class[i];
	nb_graphics_color c = colors[id_class % *N_colors];

	nb_graphics_set_source(g, c);
}

void nb_mshpoly_fill_nodes(const void *msh,
			   nb_graphics_context_t *g)
{
	fill_nodes(msh, g, NULL, set_source_null);
}

static void fill_nodes(const nb_mshpoly_t *msh,
		       nb_graphics_context_t *g,
		       const void *source_data,
		       void (*set_source)(const nb_mshpoly_t *msh,
					  nb_graphics_context_t *g,
					  uint32_t i, void *data))
{
	uint32_t N_nodes = nb_mshpoly_get_N_nodes(msh);
	for (uint32_t i = 0; i < N_nodes; i++) {
		double x = nb_msh3trg_get_x_node(msh, i);
		double y = nb_msh3trg_get_y_node(msh, i);

		nb_graphics_set_point(g, x, y, 5);

		set_source(msh, g, i, source_data);

		nb_graphics_fill(g);
	}
}

void nb_mshpoly_fill_nodes_classes(const void *msh,
				   nb_graphics_context_t *g,
				   const uint8_t *class, uint8_t N_colors,
				   const nb_graphics_color_t *colors)
{
	void *data[3];
	data[0] = (void*)class;
	data[1] = (void*)colors;
	data[2] = &N_colors;
	fill_nodes(msh, g, data, set_source_classes);
}
