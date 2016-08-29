#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>

#include "nb/math_bot.h"
#include "nb/container_bot.h"
#include "nb/graphics_bot.h"

#include "nb/geometric_bot/mesh/partition/info.h"
#include "nb/geometric_bot/mesh/partition/elements2D/msh3trg.h"

typedef struct {
	const void *msh3trg;
	uint32_t k_part;
	const uint32_t *part;
	uint32_t k_to_draw;
	double scale;
} part_data_t;

static void draw_sgm(const nb_msh3trg_t *msh3trg,
		     nb_graphics_context_t *g,
		     uint32_t isgm);

static void fill_elems(const nb_msh3trg_t *msh,
		       nb_graphics_context_t *g,
		       const void *source_data,
		       void (*set_source)(const nb_msh3trg_t *msh,
					  nb_graphics_context_t *g,
					  uint32_t i, void *data));
static void set_source_null(const nb_msh3trg_t *msh,
			    nb_graphics_context_t *g,
			    uint32_t i, void *data);

static void set_source_field_on_nodes(const nb_msh3trg_t *msh,
				      nb_graphics_context_t *g,
				      uint32_t i, void *data);

static void set_source_field_on_elems(const nb_msh3trg_t *msh,
				      nb_graphics_context_t *g,
				      uint32_t i, void *data);

static void set_source_classes(const nb_msh3trg_t *msh,
			       nb_graphics_context_t *g,
			       uint32_t i, void *data);

static void fill_nodes(const nb_msh3trg_t *msh,
		       nb_graphics_context_t *g,
		       const void *source_data,
		       void (*set_source)(const nb_msh3trg_t *msh,
					  nb_graphics_context_t *g,
					  uint32_t i, void *data));

static void calculate_partition_centers(const void *const msh3trg,
					const uint32_t *const part,
					uint32_t kpart, double **pcenter);
static void draw_triangle_partition_edge(nb_graphics_context_t *g,
					 const void *const msh3trg,
					 const uint32_t *const part,
					 uint32_t kpart, double **pcenter,
					 double scale_partitions,
					 int trg_id, int edge_id);
static void draw_triangle_partition_border(nb_graphics_context_t *g,
					   const void *const msh3trg,
					   const uint32_t *const part,
					   uint32_t kpart, double **pcenter,
					   double scale_partitions,
					   uint32_t trg_id);

static void draw_msh3trg_partition(nb_graphics_context_t *g,
				   int width, int height,
				   const void *const part_data);

static void scale_vtx(double vtx[2], double center[2], double zoom);


void nb_msh3trg_draw_wires(const void *msh,
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

void nb_msh3trg_draw_boundaries(const void *msh,
				nb_graphics_context_t *g)
{
	uint32_t N_sgm = nb_msh3trg_get_N_insgm(msh);
	for (uint32_t i = 0; i < N_sgm; i++)
		draw_sgm(msh, g, i);
}

static void draw_sgm(const nb_msh3trg_t *msh3trg,
		     nb_graphics_context_t *g,
		     uint32_t isgm)
{
	uint32_t N_vtx = nb_msh3trg_insgm_get_N_nodes(msh3trg, isgm);
	if (0 < N_vtx) {
		uint32_t nid = nb_msh3trg_insgm_get_node(msh3trg, isgm, 0);
		double x = nb_msh3trg_node_get_x(msh3trg, nid);
		double y = nb_msh3trg_node_get_y(msh3trg, nid);
		nb_graphics_move_to(g, x, y);

		for (uint32_t j = 1; j < N_vtx; j++) {
			nid = nb_msh3trg_insgm_get_node(msh3trg, isgm, j);
			x = nb_msh3trg_node_get_x(msh3trg, nid);
			y = nb_msh3trg_node_get_y(msh3trg, nid);
			nb_graphics_line_to(g, x, y);
		}
		nb_graphics_stroke(g);
	}
}

void nb_msh3trg_fill_elems(const void *msh,
			   nb_graphics_context_t *g,
			   const void* source_data,
			   void (*set_source)(const void *msh,
					      void *g,
					      uint32_t elem_id,
					      const void *data))
{
	fill_elems(msh, g, NULL, set_source_null);
}

static void set_source_null(const nb_msh3trg_t *msh,
			    nb_graphics_context_t *g,
			    uint32_t i, void *data)
{
	;/* NULL statement */
}

static void fill_elems(const nb_msh3trg_t *msh,
		       nb_graphics_context_t *g,
		       const void *source_data,
		       void (*set_source)(const nb_msh3trg_t *msh,
					  nb_graphics_context_t *g,
					  uint32_t i, void *data))
{
	uint32_t N_elems = nb_msh3trg_get_N_elems(msh);
	for (uint32_t i = 0; i < N_elems; i++) {
		uint32_t nid = nb_msh3trg_elem_get_adj(msh, i, 0);
		double x = nb_msh3trg_get_x_node(msh, nid);
		double y = nb_msh3trg_get_y_node(msh, nid);
		nb_graphics_move_to(g, x, y);

		nid = nb_msh3trg_elem_get_adj(msh, i, 1);
		x = nb_msh3trg_get_x_node(msh, nid);
		y = nb_msh3trg_get_y_node(msh, nid);
		nb_graphics_line_to(g, x, y);

		nid = nb_msh3trg_elem_get_adj(msh, i, 2);
		x = nb_msh3trg_get_x_node(msh, nid);
		y = nb_msh3trg_get_y_node(msh, nid);
		nb_graphics_line_to(g, x, y);

		nb_graphics_close_path(g);

		set_source(msh, g, i, source_data);

		nb_graphics_fill(g);
	}
}


void nb_msh3trg_fill_elems_field_on_nodes(const void *msh,
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

static void set_source_field_on_nodes(const nb_msh3trg_t *msh,
				      nb_graphics_context_t *g,
				      uint32_t i, void *data)
{
	void **cls_data = data;
	double *field = data[0];
	nb_graphics_palette_t *palette = data[1];

	uint32_t n1 = nb_msh3trg_elem_get_adj(msh, i, 0);
	uint32_t n2 = nb_msh3trg_elem_get_adj(msh, i, 1);
	uint32_t n3 = nb_msh3trg_elem_get_adj(msh, i, 2);

	double x1 = nb_msh3trg_node_get_x(msh, n1);
	double y1 = nb_msh3trg_node_get_y(msh, n1);

	double x2 = nb_msh3trg_node_get_x(msh, n2);
	double y2 = nb_msh3trg_node_get_y(msh, n2);

	double x3 = nb_msh3trg_node_get_x(msh, n3);
	double y3 = nb_msh3trg_node_get_y(msh, n3); 

	uint8_t c1[4], c2[4], c3[4];
	nb_graphics_palette_get_rgba(palette, field[n1], c1);
	nb_graphics_palette_get_rgba(palette, field[n2], c2);
	nb_graphics_palette_get_rgba(palette, field[n3], c3);

	nb_graphics_set_source_trg(g, x1, y1, x2, y2,
				   x3, y3, c1, c2, c3);
}

void nb_msh3trg_fill_elems_field_on_elems(const void *msh,
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

static void set_source_field_on_elems(const nb_msh3trg_t *msh,
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

void nb_msh3trg_fill_elems_classes(const void *msh,
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

static void set_source_classes(const nb_msh3trg_t *msh,
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

void nb_msh3trg_fill_nodes(const void *msh,
			   nb_graphics_context_t *g)
{
	fill_nodes(msh, g, NULL, set_source_null);
}

static void fill_nodes(const nb_msh3trg_t *msh,
		       nb_graphics_context_t *g,
		       const void *source_data,
		       void (*set_source)(const nb_msh3trg_t *msh,
					  nb_graphics_context_t *g,
					  uint32_t i, void *data))
{
	uint32_t N_nodes = nb_msh3trg_get_N_nodes(msh);
	for (uint32_t i = 0; i < N_nodes; i++) {
		double x = nb_msh3trg_node_get_x(msh, i);
		double y = nb_msh3trg_node_get_y(msh, i);

		nb_graphics_set_point(g, x, y, 5);

		set_source(msh, g, i, source_data);

		nb_graphics_fill(g);
	}
}

void nb_msh3trg_fill_nodes_classes(const void *msh,
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

static void calculate_partition_centers
                     (const void *const restrict msh3trg,
		      const uint32_t *const restrict part,
		      uint32_t kpart, double **pcenter)
{
	for (uint32_t k = 0; k < kpart; k++)
		memset(pcenter[k], 0, 2 * sizeof(double));

	uint32_t* ksize = calloc(kpart, sizeof(*ksize));
	uint32_t N_elems = nb_msh3trg_get_N_elems(msh3trg);
	for (uint32_t i = 0; i < N_elems; i++) {
		double xc = 0.0;
		double yc = 0.0;
		for (uint32_t j = 0; j < 3; j++) {
			uint32_t n = nb_msh3trg_elem_get_adj(msh3trg, i, j);;
			xc += nb_msh3trg_get_x_node(msh3trg, n);
			yc += nb_msh3trg_get_y_node(msh3trg, n);
		}
		pcenter[part[i]][0] += xc / 3.0;
		pcenter[part[i]][1] += yc / 3.0;
		ksize[part[i]] += 1;
	}

	for (uint32_t k = 0; k < kpart; k++) {
		pcenter[k][0] /= ksize[k];
		pcenter[k][1] /= ksize[k];
	}
}

static void draw_triangle_partition_edge(nb_graphics_context_t *g,
					 const void *const msh3trg,
					 const uint32_t *const part,
					 uint32_t kpart, double **pcenter,
					 double scale_partitions,
					 int trg_id, int edge_id)
{
	double vtx[2];
	uint32_t n1 = nb_msh3trg_elem_get_adj(msh3trg, trg_id, edge_id);
	uint32_t n2 = nb_msh3trg_elem_get_adj(msh3trg, trg_id, (edge_id+1)%3);
	vtx[0] = nb_msh3trg_get_x_node(msh3trg, n1);
	vtx[1] = nb_msh3trg_get_y_node(msh3trg, n1);

	if (scale_partitions > 0.0 && scale_partitions < 1.0)
		scale_vtx(vtx, pcenter[part[trg_id]], scale_partitions);

	nb_graphics_move_to(g, vtx[0], vtx[1]);

	vtx[0] = nb_msh3trg_get_x_node(msh3trg, n2);
	vtx[1] = nb_msh3trg_get_y_node(msh3trg, n2);

	if (scale_partitions > 0.0 && scale_partitions < 1.0)
		scale_vtx(vtx, pcenter[part[trg_id]], scale_partitions);

	nb_graphics_line_to(g, vtx[0], vtx[1]);
	nb_graphics_stroke(g);
}

static void draw_triangle_partition_border(nb_graphics_context_t *g,
					   const void *const msh3trg,
					   const uint32_t *const part,
					   uint32_t kpart, double **pcenter,
					   double scale_partitions,
					   uint32_t trg_id)
{
	uint32_t n1 = nb_msh3trg_elem_get_ngb(msh3trg, trg_id, 0);
	uint32_t n2 = nb_msh3trg_elem_get_ngb(msh3trg, trg_id, 1);
	uint32_t n3 = nb_msh3trg_elem_get_ngb(msh3trg, trg_id, 2);

	if (part[trg_id] != part[n1])
		draw_triangle_partition_edge(g, msh3trg, part, kpart, pcenter,
					     scale_partitions, trg_id, 0);

	if (part[trg_id] != part[n2])
		draw_triangle_partition_edge(g, msh3trg, part, kpart, pcenter,
					     scale_partitions, trg_id, 1);

	if (part[trg_id] != part[n3])
		draw_triangle_partition_edge(g, msh3trg, part, kpart, pcenter,
					     scale_partitions, trg_id, 2);
}

static inline void scale_vtx(double vtx[2], double center[2], double zoom)
{
	vtx[0] = (vtx[0] - center[0]) * zoom + center[0];
	vtx[1] = (vtx[1] - center[1]) * zoom + center[1];
}

void nb_msh3trg_draw_subdomain(const void *const msh3trg,
				 const char* filename, int width, int height,
				 uint32_t k_part, const uint32_t *const part,
				 uint32_t k_to_draw, double scale_partitions)
{
	if (k_part < 2) {
		nb_graphics_export(filename, width, height,
				   draw_msh3trg, msh3trg);
	} else {
		part_data_t part_data;
		part_data.msh3trg = msh3trg;
		part_data.k_part = k_part;
		part_data.part = part;
		part_data.k_to_draw = k_to_draw;
		part_data.scale = scale_partitions;
		nb_graphics_export(filename, width, height,
				   draw_msh3trg_partition,
				   &part_data);
	}
}

static void draw_msh3trg_partition(nb_graphics_context_t *g,
				   int width, int height,
				   const void *const part_data)
{
	const part_data_t *const data = part_data;
	const void *const msh3trg = data->msh3trg;

	/* Compute cam.center and cam.zoom */
	double box[4];
	nb_msh3trg_get_enveloping_box(msh3trg, box);

	nb_graphics_enable_camera(g);
	nb_graphics_camera_t* cam = nb_graphics_get_camera(g);
	nb_graphics_cam_fit_box(cam, box, width, height);

	/* Calculate center of partitions */
	double** pcenter = NULL;
	if (data->scale > 0.0 && data->scale < 1.0) {
		pcenter = malloc(data->k_part * sizeof(*pcenter));
		for (uint32_t k = 0; k < data->k_part; k++)
			pcenter[k] = calloc(2, sizeof(*(pcenter[k])));
		calculate_partition_centers(msh3trg, data->part,
					    data->k_part, pcenter);
	}

	/* Draw triangles */
	uint32_t N_elems = nb_msh3trg_get_N_elems(msh3trg);
	for (uint32_t i = 0; i < N_elems; i++) {
		if (data->part[i] != data->k_to_draw &&
		    data->k_to_draw < data->k_part)
			continue;
    
		uint32_t n1 = nb_msh3trg_elem_get_adj(msh3trg, i, 0);
		uint32_t n2 = nb_msh3trg_elem_get_adj(msh3trg, i, 1);
		uint32_t n3 = nb_msh3trg_elem_get_adj(msh3trg, i, 2);

		double vtx[2];
		vtx[0] = nb_msh3trg_get_x_node(msh3trg, n1);
		vtx[1] = nb_msh3trg_get_y_node(msh3trg, n1);

		if (data->scale > 0.0 && data->scale < 1.0)
			scale_vtx(vtx, pcenter[data->part[i]], data->scale);

		nb_graphics_move_to(g, vtx[0], vtx[1]);

		vtx[0] = nb_msh3trg_get_x_node(msh3trg, n2);
		vtx[1] = nb_msh3trg_get_y_node(msh3trg, n2);
		if (data->scale > 0.0 && data->scale < 1.0)
			scale_vtx(vtx, pcenter[data->part[i]], data->scale);
		nb_graphics_line_to(g, vtx[0], vtx[1]);

		vtx[0] = nb_msh3trg_get_x_node(msh3trg, n3);
		vtx[1] = nb_msh3trg_get_y_node(msh3trg, n3);
		if (data->scale > 0.0 && data->scale < 1.0)
			scale_vtx(vtx, pcenter[data->part[i]], data->scale);

		nb_graphics_line_to(g, vtx[0], vtx[1]);
		nb_graphics_close_path(g);

		nb_graphics_set_source(g, NB_LIGHT_GRAY);
		nb_graphics_fill_preserve(g);

		nb_graphics_set_source_rgb(g, 75, 128, 255);
		nb_graphics_set_line_width(g, 0.5);

		nb_graphics_stroke(g);

		/* Draw border */
		nb_graphics_set_source_rgb(g, 255, 0, 200);
		nb_graphics_set_line_width(g, 1.0);

		draw_triangle_partition_border(g, msh3trg, data->part,
					       data->k_part,
					       pcenter, data->scale, i);
	}

	if (data->scale > 0.0 && data->scale < 1.0) {
		for (uint32_t k = 0; k < data->k_part; k++)
			free(pcenter[k]);
		free(pcenter);
	}

	/* Draw input segments */
	if (data->k_to_draw >= data->k_part) {
		nb_graphics_set_source_rgb(g, 255, 0, 200);

		nb_graphics_set_line_width(g, 1.0);
		uint32_t N_sgm = nb_msh3trg_get_N_insgm(msh3trg);
		for (uint32_t i = 0; i < N_sgm; i++) {
			uint32_t N_vtx = nb_msh3trg_get_N_nodes_x_insgm(msh3trg,
									i);
			if (0 < N_vtx) {
				uint32_t n1 = nb_msh3trg_get_node_x_insgm(msh3trg, i, 0);
				nb_graphics_move_to(g,
						    nb_msh3trg_get_x_node(msh3trg, n1),
						    nb_msh3trg_get_y_node(msh3trg, n1));
				for (uint32_t j = 0; j < N_vtx; j++) {
					uint32_t n2 = nb_msh3trg_get_node_x_insgm(msh3trg, i, j);
					nb_graphics_line_to(g,
							    nb_msh3trg_get_x_node(msh3trg, n2),
							    nb_msh3trg_get_y_node(msh3trg, n2));
				}
				nb_graphics_stroke(g);
			}
		}
	}
}
