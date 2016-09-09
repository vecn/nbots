#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>

#include "nb/graphics_bot.h"

#include "nb/memory_bot.h"
#include "nb/container_bot.h"
#include "nb/geometric_bot/utils2D.h"
#include "nb/geometric_bot/mesh/partition.h"
#include "nb/geometric_bot/mesh/partition/info.h"
#include "nb/geometric_bot/mesh/partition/draw.h"
#include "partition_struct.h"

#define COLOR_EDGE NB_DARK_GRAY
#define COLOR_ELEM NB_LIGHT_BLUE
#define COLOR_SGM NB_BLACK
#define PALETTE_FIELD NB_RAINBOW
#define COLOR_CLASS_0 NB_BLUE
#define COLOR_CLASS_1 NB_RED
#define COLOR_CLASS_2 NB_VIOLET
#define COLOR_CLASS_3 NB_AZURE
#define COLOR_CLASS_4 NB_GREEN
#define COLOR_CLASS_5 NB_YELLOW
#define COLOR_CLASS_6 NB_ROSE
#define COLOR_CLASS_7 NB_CHARTREUSE
#define COLOR_CLASS_8 NB_LIGHT_GRAY
#define COLOR_CLASS_9 NB_AQUAMARIN

#define PALETTE_W 30
#define PALETTE_H 250
#define PALETTE_MARGIN 40

typedef struct {
	const nb_partition_t *part;
	const void *values;
	nb_partition_entity vals_entity;
	nb_partition_array_type vals_type;
	bool draw_wires;
} draw_data;

static void init_draw_data(draw_data *data,
			   const nb_partition_t *part,
			   nb_partition_entity vals_entity,
			   nb_partition_array_type vals_type,
			   const void *values,
			   bool draw_wires);

static void draw(nb_graphics_context_t *g, int width, int height,
		 const void *draw_data);
static void set_camera(nb_graphics_context_t *g, int width, int height,
		       const nb_partition_t *part);
static void fill(const nb_partition_t *part,
		 nb_graphics_context_t *g,
		 const draw_data *data);

static void fill_elems_field_on_nodes(const nb_partition_t *part,
				      nb_graphics_context_t *g,
				      const double *values);
static void normalize_values(double *normalized_values,
			     const double *values, uint32_t N);
static void fill_elems_field_on_elems(const nb_partition_t *part,
				      nb_graphics_context_t *g,
				      const double *values);
static void fill_elems_classes(const nb_partition_t *part,
			       nb_graphics_context_t *g,
			       const uint8_t *class);
static void set_class_colors(nb_graphics_color_t color[10]);
static void fill_elems(const nb_partition_t *part,
		       nb_graphics_context_t *g);

static void draw_wires(const nb_partition_t *part,
		       nb_graphics_context_t *g);

static void draw_boundaries(const nb_partition_t *part,
			    nb_graphics_context_t *g);

static void fill_nodes_classes(const nb_partition_t *part,
			       nb_graphics_context_t *g,
			       const uint8_t *class);
static void add_palette(const nb_partition_t *part,
			const draw_data *data,
			nb_graphics_context_t *g, int width, int height);
static void get_min_max(const nb_partition_t *part,
			const draw_data *data,
			double *min, double *max);

void nb_partition_export_draw(const nb_partition_t *part,
			      const char *filename,
			      int width, int height,
			      nb_partition_entity vals_entity,
			      nb_partition_array_type vals_type,
			      const void *values,
			      bool draw_wires)
{
	draw_data data;
	init_draw_data(&data, part, vals_entity,
		       vals_type, values, draw_wires);

	nb_graphics_export(filename, width, height,
			   draw, &data);
}

static void init_draw_data(draw_data *data,
			   const nb_partition_t *part,
			   nb_partition_entity vals_entity,
			   nb_partition_array_type vals_type,
			   const void *values,
			   bool draw_wires)
{
	data->part = part;
	data->values = values;
	data->vals_entity = vals_entity;
	data->vals_type = vals_type;
	data->draw_wires = draw_wires;
}

static void draw(nb_graphics_context_t *g, int width, int height,
		 const void *data_ptr)
{
	const draw_data *data = data_ptr;
	const nb_partition_t *part = data->part;

	if (!nb_graphics_is_camera_enabled(g))
		set_camera(g, width, height, part);

	fill(part, g, data);

	if (data->draw_wires)
		draw_wires(part, g);

	draw_boundaries(part, g);

	if (NB_NODE == data->vals_entity && NB_CLASS == data->vals_type)
		fill_nodes_classes(part, g, data->values);

	if (NB_FIELD == data->vals_type)
		add_palette(part, data, g, width, height);
}

static void set_camera(nb_graphics_context_t *g, int width, int height,
		       const nb_partition_t *part)
{
	double box[4];
	part->get_enveloping_box(part->msh, box);

	nb_graphics_enable_camera(g);
	nb_graphics_camera_t* cam = nb_graphics_get_camera(g);
	nb_graphics_cam_fit_box(cam, box, width, height);
}

static void fill(const nb_partition_t *part,
		 nb_graphics_context_t *g,
		 const draw_data *data)
{
	nb_partition_entity enty = data->vals_entity;
	nb_partition_array_type type = data->vals_type;

	if (NB_NODE == enty && NB_FIELD == type)
		fill_elems_field_on_nodes(part, g, data->values);

	else if (NB_ELEMENT == enty && NB_FIELD == type)
		fill_elems_field_on_elems(part, g, data->values);

	else if (NB_ELEMENT == enty && NB_CLASS == type)
		fill_elems_classes(part, g, data->values);

	else
		fill_elems(part, g);
}

static void fill_elems_field_on_nodes(const nb_partition_t *part,
				      nb_graphics_context_t *g,
				      const double *values)
{
	uint32_t N = part->get_N_nodes(part->msh);
	uint32_t memsize = N * sizeof(double);
	double *normalized_values = NB_SOFT_MALLOC(memsize);

	normalize_values(normalized_values, values, N);
	
	part->graphics.fill_elems_field_on_nodes(part->msh, g,
						 normalized_values,
						 PALETTE_FIELD);

	NB_SOFT_FREE(memsize, normalized_values);
}

static void normalize_values(double *normalized_values,
			     const double *values, uint32_t N)
{
	uint32_t min_id;
	uint32_t max_id;
	vcn_array_get_min_max_ids(values, N, sizeof(*values),
				  vcn_compare_double,
				  &(min_id), &(max_id));

	double min = values[min_id];
	double range = values[max_id] - min;
	for (uint32_t i = 0; i < N; i++)
		normalized_values[i] = (values[i] - min) / range;
}

static void fill_elems_field_on_elems(const nb_partition_t *part,
				      nb_graphics_context_t *g,
				      const double *values)
{
	uint32_t N = part->get_N_elems(part->msh);
	uint32_t memsize = N * sizeof(double);
	double *normalized_values = NB_SOFT_MALLOC(memsize);

	normalize_values(normalized_values, values, N);
	
	part->graphics.fill_elems_field_on_elems(part->msh, g,
						 normalized_values,
						 PALETTE_FIELD);

	NB_SOFT_FREE(memsize, normalized_values);
}

static void fill_elems_classes(const nb_partition_t *part,
			       nb_graphics_context_t *g,
			       const uint8_t *class)
{
	nb_graphics_color_t color[10];
	set_class_colors(color);
	part->graphics.fill_elems_classes(part->msh, g, class, 10, color);
}

static void set_class_colors(nb_graphics_color_t color[10])
{
	color[0] = COLOR_CLASS_0;
	color[1] = COLOR_CLASS_1;
	color[2] = COLOR_CLASS_2;
	color[3] = COLOR_CLASS_3;
	color[4] = COLOR_CLASS_4;
	color[5] = COLOR_CLASS_5;
	color[6] = COLOR_CLASS_6;
	color[7] = COLOR_CLASS_7;
	color[8] = COLOR_CLASS_8;
	color[9] = COLOR_CLASS_9;
}

static void fill_elems(const nb_partition_t *part,
		       nb_graphics_context_t *g)
{
	nb_graphics_set_source(g, COLOR_ELEM);
	part->graphics.fill_elems(part->msh, g);
}

static void draw_wires(const nb_partition_t *part,
		       nb_graphics_context_t *g)
{
	nb_graphics_set_source(g, COLOR_EDGE);
	nb_graphics_set_line_width(g, 1.0);
	part->graphics.draw_wires(part->msh, g);
}

static void draw_boundaries(const nb_partition_t *part,
			    nb_graphics_context_t *g)
{
	nb_graphics_set_source(g, COLOR_SGM);
	nb_graphics_set_line_width(g, 1.5);
	part->graphics.draw_boundaries(part->msh, g);
}

static void fill_nodes_classes(const nb_partition_t *part,
			       nb_graphics_context_t *g,
			       const uint8_t *class)
{
	nb_graphics_color_t color[10];
	set_class_colors(color);
	part->graphics.fill_nodes_classes(part->msh, g, class, 10, color);
}

static void add_palette(const nb_partition_t *part,
			const draw_data *data,
			nb_graphics_context_t *g, int width, int height)
{
	nb_graphics_palette_t *palette =
		nb_graphics_palette_create_preset(PALETTE_FIELD);
	double min, max;
	get_min_max(part, data, &min, &max);
	
	float label_width = 55;
	nb_graphics_palette_draw(g, palette,
				 width - PALETTE_W - PALETTE_MARGIN -
				 label_width,
				 height - PALETTE_H - PALETTE_MARGIN,
				 PALETTE_W,
				 PALETTE_H,
				 1, min, max);
	nb_graphics_palette_destroy(palette);
}

static void get_min_max(const nb_partition_t *part,
			const draw_data *data,
			double *min, double *max)
{
	uint32_t N;
	if (NB_ELEMENT == data->vals_entity)
		N = nb_partition_get_N_elems(part);
	else
		N = nb_partition_get_N_nodes(part);

	const double *values = data->values;
	*min = values[0];
	*max = values[0];
	for (uint32_t i = 1; i < N; i++) {
		if (values[i] < *min)
			*min = values[i];
		else if (values[i] > *max)
			*max = values[i];
	}
}
