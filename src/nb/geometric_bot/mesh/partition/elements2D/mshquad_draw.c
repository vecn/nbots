#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>

#include "nb/math_bot.h"
#include "nb/container_bot.h"
#include "nb/graphics_bot.h"

#include "nb/geometric_bot/mesh/partition/elements2D/mshquad.h"
#include "nb/geometric_bot/mesh/partition/elements2D/mshquad_draw.h"

static void draw_sgm(const nb_mshquad_t *msh,
		     nb_graphics_context_t *g,
		     uint32_t isgm);
static void set_source_null(const nb_mshquad_t *msh,
			    nb_graphics_context_t *g,
			    uint32_t i, void *data);
static void fill_elems(const nb_mshquad_t *msh,
		       nb_graphics_context_t *g,
		       void *source_data,
		       void (*set_source)(const nb_mshquad_t *msh,
					  nb_graphics_context_t *g,
					  uint32_t i, void *data));
static void set_source_field_on_elems(const nb_mshquad_t *msh,
				      nb_graphics_context_t *g,
				      uint32_t i, void *data);
static void set_source_classes(const nb_mshquad_t *msh,
			       nb_graphics_context_t *g,
			       uint32_t i, void *data);

static void fill_nodes(const nb_mshquad_t *msh,
		       nb_graphics_context_t *g,
		       void *source_data,
		       void (*set_source)(const nb_mshquad_t *msh,
					  nb_graphics_context_t *g,
					  uint32_t i, void *data));

void nb_mshquad_draw_wires(const void *msh,
			   nb_graphics_context_t *g)
{
	uint32_t N_edg = nb_mshquad_get_N_edges(msh);
	for (uint32_t i = 0; i < N_edg; i++) {
		uint32_t n1 = nb_mshquad_edge_get_1n(msh, i);
		uint32_t n2 = nb_mshquad_edge_get_2n(msh, i);

		double x = nb_mshquad_node_get_x(msh, n1);
		double y = nb_mshquad_node_get_y(msh, n1);
		nb_graphics_move_to(g, x, y);

		x = nb_mshquad_node_get_x(msh, n2);
		y = nb_mshquad_node_get_y(msh, n2);
		nb_graphics_line_to(g, x, y);

		nb_graphics_stroke(g);
	}
}

void nb_mshquad_draw_boundaries(const void *msh,
				nb_graphics_context_t *g)
{
	uint32_t N_sgm = nb_mshquad_get_N_insgm(msh);
	for (uint32_t i = 0; i < N_sgm; i++)
		draw_sgm(msh, g, i);
}

static void draw_sgm(const nb_mshquad_t *msh,
		     nb_graphics_context_t *g,
		     uint32_t isgm)
{
	uint32_t N_vtx = nb_mshquad_insgm_get_N_nodes(msh, isgm);
	if (0 < N_vtx) {
		uint32_t nid = nb_mshquad_insgm_get_node(msh, isgm, 0);
		double x = nb_mshquad_node_get_x(msh, nid);
		double y = nb_mshquad_node_get_y(msh, nid);
		nb_graphics_move_to(g, x, y);

		for (uint32_t j = 1; j < N_vtx; j++) {
			nid = nb_mshquad_insgm_get_node(msh, isgm, j);
			x = nb_mshquad_node_get_x(msh, nid);
			y = nb_mshquad_node_get_y(msh, nid);
			nb_graphics_line_to(g, x, y);
		}
		nb_graphics_stroke(g);
	}
}

void nb_mshquad_fill_elems(const void *msh,
			   nb_graphics_context_t *g)
{
	fill_elems(msh, g, NULL, set_source_null);
}

static void set_source_null(const nb_mshquad_t *msh,
			    nb_graphics_context_t *g,
			    uint32_t i, void *data)
{
	;/* NULL statement */
}

static void fill_elems(const nb_mshquad_t *msh,
		       nb_graphics_context_t *g,
		       void *source_data,
		       void (*set_source)(const nb_mshquad_t *msh,
					  nb_graphics_context_t *g,
					  uint32_t i, void *data))
{
	uint32_t N_elems = nb_mshquad_get_N_elems(msh);
	for (uint32_t i = 0; i < N_elems; i++) {
		uint32_t n1 = nb_mshquad_elem_get_adj(msh, i, 0);
		uint32_t n2 = nb_mshquad_elem_get_adj(msh, i, 1);
		uint32_t n3 = nb_mshquad_elem_get_adj(msh, i, 2);
		uint32_t n4 = nb_mshquad_elem_get_adj(msh, i, 3);

		double x = nb_mshquad_node_get_x(msh, n1);
		double y = nb_mshquad_node_get_y(msh, n1);
		nb_graphics_move_to(g, x, y);

		x = nb_mshquad_node_get_x(msh, n2);
		y = nb_mshquad_node_get_y(msh, n2);
		nb_graphics_line_to(g, x, y);

		x = nb_mshquad_node_get_x(msh, n3);
		y = nb_mshquad_node_get_y(msh, n3);
		nb_graphics_line_to(g, x, y);

		if (n4 < nb_mshquad_get_N_nodes(msh)) {
			x = nb_mshquad_node_get_x(msh, n4);
			y = nb_mshquad_node_get_y(msh, n4);
			nb_graphics_line_to(g, x, y);
		}
		nb_graphics_close_path(g);

		set_source(msh, g, i, source_data);

		nb_graphics_fill(g);
	}
}

void nb_mshquad_fill_elems_field_on_nodes(const void *msh,
					  nb_graphics_context_t *g,
					  const double *normalized_field,
					  nb_graphics_palette_preset palette)
{
	nb_graphics_palette_t *pal = nb_graphics_palette_create_preset(palette);

	uint32_t N_elems = nb_mshquad_get_N_elems(msh);
	for (uint32_t i = 0; i < N_elems; i++) {
		uint32_t n1 = nb_mshquad_elem_get_adj(msh, i, 0);
		double x1 = nb_mshquad_node_get_x(msh, n1);
		double y1 = nb_mshquad_node_get_y(msh, n1);
		
		uint32_t n2 = nb_mshquad_elem_get_adj(msh, i, 1);
		double x2 = nb_mshquad_node_get_x(msh, n2);
		double y2 = nb_mshquad_node_get_y(msh, n2);

		uint32_t n3 = nb_mshquad_elem_get_adj(msh, i, 2);
		double x3 = nb_mshquad_node_get_x(msh, n3);
		double y3 = nb_mshquad_node_get_y(msh, n3);


		uint8_t c1[4], c2[4], c3[4];
		nb_graphics_palette_get_rgba(pal, normalized_field[n1], c1);
		nb_graphics_palette_get_rgba(pal, normalized_field[n2], c2);
		nb_graphics_palette_get_rgba(pal, normalized_field[n3], c3);

		nb_graphics_set_source_trg(g, x1, y1, x2, y2,
					   x3, y3, c1, c2, c3);

		nb_graphics_move_to(g, x1, y1);
		nb_graphics_line_to(g, x2, y2);
		nb_graphics_line_to(g, x3, y3);
		nb_graphics_close_path(g);

		nb_graphics_fill(g);

		uint32_t n4 = nb_mshquad_elem_get_adj(msh, i, 3);
		if (n4 < nb_mshquad_get_N_nodes(msh)) {
			double x4 = nb_mshquad_node_get_x(msh, n4);
			double y4 = nb_mshquad_node_get_y(msh, n4);

			uint8_t c4[4];
			double field4 = normalized_field[n4];
			nb_graphics_palette_get_rgba(pal, field4, c4);

			nb_graphics_set_source_trg(g, x1, y1, x3, y3,
						   x4, y4, c1, c3, c4);

			nb_graphics_move_to(g, x1, y1);
			nb_graphics_line_to(g, x3, y3);
			nb_graphics_line_to(g, x4, y4);
			nb_graphics_close_path(g);

			nb_graphics_fill(g);
		}
	}	

	nb_graphics_palette_destroy(pal);
}

void nb_mshquad_fill_elems_field_on_elems(const void *msh,
					  nb_graphics_context_t *g,
					  const double *normalized_field,
					  nb_graphics_palette_preset palette)
{
	void *data[2];
	data[0] = (void*) normalized_field;
	data[1] = nb_graphics_palette_create_preset(palette);
	
	fill_elems(msh, g, data, set_source_field_on_elems);
	
	nb_graphics_palette_destroy(data[1]);
}

static void set_source_field_on_elems(const nb_mshquad_t *msh,
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

void nb_mshquad_fill_elems_classes(const void *msh,
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

static void set_source_classes(const nb_mshquad_t *msh,
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

void nb_mshquad_fill_nodes(const void *msh,
			   nb_graphics_context_t *g)
{
	fill_nodes(msh, g, NULL, set_source_null);
}

static void fill_nodes(const nb_mshquad_t *msh,
		       nb_graphics_context_t *g,
		       void *source_data,
		       void (*set_source)(const nb_mshquad_t *msh,
					  nb_graphics_context_t *g,
					  uint32_t i, void *data))
{
	uint32_t N_nodes = nb_mshquad_get_N_nodes(msh);
	for (uint32_t i = 0; i < N_nodes; i++) {
		double x = nb_mshquad_node_get_x(msh, i);
		double y = nb_mshquad_node_get_y(msh, i);

		nb_graphics_set_point(g, x, y, 5);

		set_source(msh, g, i, source_data);

		nb_graphics_fill(g);
	}
}

void nb_mshquad_fill_nodes_classes(const void *msh,
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

