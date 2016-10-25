#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>

#include "nb/math_bot.h"
#include "nb/memory_bot.h"
#include "nb/container_bot.h"
#include "nb/graphics_bot.h"

#include "nb/geometric_bot/mesh/mesh2D/elements2D/mshpoly.h"
#include "nb/geometric_bot/mesh/mesh2D/elements2D/mshpoly_draw.h"

static void draw_sgm(const nb_mshpoly_t *poly,
		     nb_graphics_context_t *g,
		     uint32_t isgm);
static void set_source_null(const nb_mshpoly_t *msh,
			    nb_graphics_context_t *g,
			    uint32_t i, void *data);
static void fill_elems(const nb_mshpoly_t *poly,
		       nb_graphics_context_t *g,
		       void *source_data,
		       void (*set_source)(const nb_mshpoly_t *msh,
					  nb_graphics_context_t *g,
					  uint32_t i, void *data));
static void get_centroid_color(const void *msh,
			       const nb_palette_t *pal,
			       const double *field, uint32_t elem_id,
			       uint8_t cc[4]);
static void set_source_field_on_elems(const nb_mshpoly_t *msh,
				      nb_graphics_context_t *g,
				      uint32_t i, void *data);
static void set_source_classes(const nb_mshpoly_t *msh,
			       nb_graphics_context_t *g,
			       uint32_t i, void *data);
static void fill_nodes(const nb_mshpoly_t *msh,
		       nb_graphics_context_t *g,
		       void *source_data,
		       void (*set_source)(const nb_mshpoly_t *msh,
					  nb_graphics_context_t *g,
					  uint32_t i, void *data));
static double get_centroid_val(const void *msh,
			       const double *field, uint32_t elem_id);
static void draw_trg_level_set_intersection(nb_graphics_context_t *g,
					    const double t1[2],
					    const double t2[2],
					    const double t3[2],
					    double v1, double v2, double v3,
					    double level_set);

void nb_mshpoly_draw_wires(const void *msh,
			   nb_graphics_context_t *g)
{
	uint32_t N_edg = nb_mshpoly_get_N_edges(msh);
	for (uint32_t i = 0; i < N_edg; i++) {
		uint32_t n1 = nb_mshpoly_edge_get_1n(msh, i);
		uint32_t n2 = nb_mshpoly_edge_get_2n(msh, i);

		double x = nb_mshpoly_node_get_x(msh, n1);
		double y = nb_mshpoly_node_get_y(msh, n1);
		nb_graphics_move_to(g, x, y);

		x = nb_mshpoly_node_get_x(msh, n2);
		y = nb_mshpoly_node_get_y(msh, n2);
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
		uint32_t nid = nb_mshpoly_insgm_get_node(poly, isgm, 0);
		double x = nb_mshpoly_node_get_x(poly, nid);
		double y = nb_mshpoly_node_get_y(poly, nid);
		nb_graphics_move_to(g, x, y);

		for (uint32_t j = 1; j < N_vtx; j++) {
			nid = nb_mshpoly_insgm_get_node(poly, isgm, j);
			x = nb_mshpoly_node_get_x(poly, nid);
			y = nb_mshpoly_node_get_y(poly, nid);
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
		       void *source_data,
		       void (*set_source)(const nb_mshpoly_t *msh,
					  nb_graphics_context_t *g,
					  uint32_t i, void *data))
{
	uint32_t N_elems = nb_mshpoly_get_N_elems(poly);
	for (uint32_t i = 0; i < N_elems; i++) {
		uint32_t id = nb_mshpoly_elem_get_adj(poly, i, 0);
		double x = nb_mshpoly_node_get_x(poly, id);
		double y = nb_mshpoly_node_get_y(poly, id);
		nb_graphics_move_to(g, x, y);

		uint16_t N_adj = nb_mshpoly_elem_get_N_adj(poly, i);
		for (uint16_t j = 1; j < N_adj; j++) {
			id = nb_mshpoly_elem_get_adj(poly, i, j);
			x = nb_mshpoly_node_get_x(poly, id);
			y = nb_mshpoly_node_get_y(poly, id);
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
					  nb_palette_preset palette)
{
	nb_palette_t *pal = 
		nb_palette_create_preset(palette);

	uint32_t N_elems = nb_mshpoly_get_N_elems(msh);
	for (uint32_t i = 0; i < N_elems; i++) {
		double xc = nb_mshpoly_elem_get_x(msh, i);
		double yc = nb_mshpoly_elem_get_y(msh, i);
		uint8_t cc[4];
		get_centroid_color(msh, pal, normalized_field, i, cc);

		uint16_t N_adj = nb_mshpoly_elem_get_N_adj(msh, i);
		for (uint16_t j = 0; j < N_adj; j++) {
			uint32_t n1 = nb_mshpoly_elem_get_adj(msh, i, j);
			double x1 = nb_mshpoly_node_get_x(msh, n1);
			double y1 = nb_mshpoly_node_get_y(msh, n1);

			uint32_t n2 = nb_mshpoly_elem_get_adj(msh, i,
							      (j+1) % N_adj);
			double x2 = nb_mshpoly_node_get_x(msh, n2);
			double y2 = nb_mshpoly_node_get_y(msh, n2);

			nb_graphics_move_to(g, xc, yc);
			nb_graphics_line_to(g, x1, y1);
			nb_graphics_line_to(g, x2, y2);
			nb_graphics_close_path(g);
			
			uint8_t c1[4], c2[4];
			nb_palette_get_rgba(pal, normalized_field[n1],
						     c1);
			nb_palette_get_rgba(pal, normalized_field[n2],
						     c2);

			nb_graphics_set_source_trg(g, xc, yc, x1, y1,
						   x2, y2, cc, c1, c2);
			nb_graphics_fill(g);
		}
	}
	nb_palette_destroy(pal);
}

static void get_centroid_color(const void *msh,
			       const nb_palette_t *pal,
			       const double *field, uint32_t elem_id,
			       uint8_t cc[4])
{
	uint16_t rgba_sum[4] = {0, 0, 0, 0};
	uint16_t N_adj = nb_mshpoly_elem_get_N_adj(msh, elem_id);
	for (uint16_t j = 0; j < N_adj; j++) {
		uint32_t nid = nb_mshpoly_elem_get_adj(msh, elem_id, j);
		uint8_t rgba[4];
		nb_palette_get_rgba(pal, field[nid], rgba);
			
		rgba_sum[0] += rgba[0];
		rgba_sum[1] += rgba[1];
		rgba_sum[2] += rgba[2];
		rgba_sum[3] += rgba[3];
	}
	cc[0] = rgba_sum[0] / N_adj;
	cc[1] = rgba_sum[1] / N_adj;
	cc[2] = rgba_sum[2] / N_adj;
	cc[3] = rgba_sum[3] / N_adj;
}

void nb_mshpoly_fill_elems_field_on_elems(const void *msh,
					  nb_graphics_context_t *g,
					  const double *normalized_field,
					  nb_palette_preset palette)
{
	void *data[2];
	data[0] = (void*) normalized_field;
	data[1] = nb_palette_create_preset(palette);
	
	fill_elems(msh, g, data, set_source_field_on_elems);
	
	nb_palette_destroy(data[1]);
}

static void set_source_field_on_elems(const nb_mshpoly_t *msh,
				      nb_graphics_context_t *g,
				      uint32_t i, void *data)
{
	void **cls_data = data;
	double *field = cls_data[0];
	nb_palette_t *palette = cls_data[1];
	
	uint8_t c[4];
	nb_palette_get_rgba(palette, field[i], c);
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
	nb_graphics_color_t c = colors[id_class % *N_colors];

	nb_graphics_set_source(g, c);
}

void nb_mshpoly_fill_nodes(const void *msh,
			   nb_graphics_context_t *g)
{
	fill_nodes(msh, g, NULL, set_source_null);
}

static void fill_nodes(const nb_mshpoly_t *msh,
		       nb_graphics_context_t *g,
		       void *source_data,
		       void (*set_source)(const nb_mshpoly_t *msh,
					  nb_graphics_context_t *g,
					  uint32_t i, void *data))
{
	uint32_t N_nodes = nb_mshpoly_get_N_nodes(msh);
	for (uint32_t i = 0; i < N_nodes; i++) {
		double x = nb_mshpoly_node_get_x(msh, i);
		double y = nb_mshpoly_node_get_y(msh, i);

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

void nb_mshpoly_draw_level_set(const void *msh,
			       nb_graphics_context_t *g,
			       const double *field_on_nodes,
			       double level_set)
{
	uint32_t N_elems = nb_mshpoly_get_N_elems(msh);
	for (uint32_t i = 0; i < N_elems; i++) {
		double xc[2];
		xc[0] = nb_mshpoly_elem_get_x(msh, i);
		xc[1] = nb_mshpoly_elem_get_y(msh, i);
		double vc = get_centroid_val(msh, field_on_nodes, i);

		uint16_t N_adj = nb_mshpoly_elem_get_N_adj(msh, i);
		for (uint16_t j = 0; j < N_adj; j++) {
			uint32_t n1 = nb_mshpoly_elem_get_adj(msh, i, j);
			double x1[2];
			x1[0] = nb_mshpoly_node_get_x(msh, n1);
			x1[1] = nb_mshpoly_node_get_y(msh, n1);
			double v1 = field_on_nodes[n1];

			uint32_t n2 = nb_mshpoly_elem_get_adj(msh, i,
							      (j+1) % N_adj);
			double x2[2];
			x2[0] = nb_mshpoly_node_get_x(msh, n2);
			x2[1] = nb_mshpoly_node_get_y(msh, n2);
			double v2 = field_on_nodes[n2];

			draw_trg_level_set_intersection(g, x1, x2, xc,
							v1, v2, vc, level_set);
		}
	}
}

static double get_centroid_val(const void *msh,
			       const double *field, uint32_t elem_id)
{
	uint16_t N = nb_mshpoly_elem_get_N_adj(msh, elem_id);

	uint32_t memsize = N * 3 * sizeof(double);
	char *memblock = nb_soft_allocate_mem(memsize);
	double *w = (void*) memblock;
	double *xi = (void*) (memblock + N * sizeof(double));

	for (uint16_t j = 0; j < N; j++) {
		uint32_t nid = nb_mshpoly_elem_get_adj(msh, elem_id, j);
		xi[j * 2] = nb_mshpoly_node_get_x(msh, nid);
		xi[j*2+1] = nb_mshpoly_node_get_y(msh, nid);
	}
	double xc[2];
	xc[0] = nb_mshpoly_elem_get_x(msh, elem_id);
	xc[1] = nb_mshpoly_elem_get_y(msh, elem_id);

	nb_nonpolynomial_eval(N, 2, xi, NULL, xc, w);

	double val = 0;
	for (uint16_t j = 0; j < N; j++) {
		uint32_t nid = nb_mshpoly_elem_get_adj(msh, elem_id, j);
		val += w[j] * field[nid];
	}

	nb_soft_free_mem(memsize, memblock);
	return val;
}

static void draw_trg_level_set_intersection(nb_graphics_context_t *g,
					    const double t1[2],
					    const double t2[2],
					    const double t3[2],
					    double v1, double v2, double v3,
					    double level_set)
{
	if (nb_utils2D_level_set_intersects_trg(v1, v2, v3, level_set)) {
		double a[2], b[2];
		nb_utils2D_get_trg_level_set_intersection(t1, t2, t3, v1, v2,
							  v3, level_set, a, b);
		nb_graphics_move_to(g, a[0], a[1]);
		nb_graphics_line_to(g, b[0], b[1]);
		nb_graphics_stroke(g);
	}
}
