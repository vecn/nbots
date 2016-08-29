#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>

#include "nb/graphics_bot.h"

#include "nb/memory_bot.h"
#include "nb/geometric_bot/mesh/partition.h"
#include "nb/geometric_bot/mesh/partition/info.h"
#include "nb/geometric_bot/mesh/partition/draw.h"
#include "partition_struct.h"

typedef struct {
	const nb_partition_t *part;
	const void *values;
	nb_partition_entity vals_entity;
	nb_partition_array_type vals_type;
	bool draw_wires;
} draw_data;

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
			   draw, &draw_data);
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
		 const void *draw_data)
{
	const nb_partition_draw_data *data = draw_data;
	const nb_partition_t *part = data->part;

	if (!nb_graphics_is_camera_enabled(g))
		set_camera(g, width, height, part);

	fill(part, g, data);

	if (data->draw_wires)
		draw_wires(part, g);

	draw_boundaries(part, g);

	if (NB_NODE == data->vals_entity && NB_CLASS == data->vals_type)
		fill_nodes_classes(part, g, data->values);
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
	
	part->di.fill_elems_field_on_nodes(part->msh, g,
					   normalized_values,
					   NB_PALETTE_RAINBOW);

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
	
	part->di.fill_elems_field_on_elems(part->msh, g,
					   normalized_values,
					   NB_PALETTE_RAINBOW);

	NB_SOFT_FREE(memsize, normalized_values);
}

static void fill_elems_classes(const nb_partition_t *part,
			       nb_graphics_context_t *g,
			       const uint8_t *class)
{
	nb_graphics_color_t color[10];
	set_class_colors(color);
	part->di.fill_elems_classes(part->msh, g, class, 10, color);
}

static void set_class_colors(nb_graphics_color_t color[10])
{
	color[0] = NB_BLUE;
	color[1] = NB_RED;
	color[2] = NB_VIOLET;
	color[3] = NB_AZURE;
	color[4] = NB_GREEN;
	color[5] = NB_YELLOW;
	color[6] = NB_ROSE;
	color[7] = NB_CHARTREUSE;
	color[8] = NB_LIGHT_GRAY;
	color[9] = NB_AQUAMARIN;
}

static void fill_elems(const nb_partition_t *part,
		       nb_graphics_context_t *g)
{
	nb_graphics_set_color(g, NB_LIGHT_BLUE);
	part->di.fill_elems(part->msh, g);
}

static void draw_wires(const nb_partition_t *part,
		       nb_graphics_context_t *g)
{
	nb_graphics_set_source(g, NB_LIGHT_PURPLE);
	nb_graphics_set_line_width(g, 0.5);
	part->di.draw_wires(part->msh, g);
}

static void draw_boundaries(const nb_partition_t *part,
			    nb_graphics_context_t *g)
{
	nb_graphics_set_source(g, NB_PURPLE);
	nb_graphics_set_line_width(g, 1.5);
	part->di.draw_boundaries(part->msh, g);
}

static void fill_nodes_classes(const nb_partition_t *part,
			       nb_graphics_context_t *g,
			       const uint8_t *class)
{
	nb_graphics_color_t color[10];
	set_class_colors(color);
	part->di.fill_nodes_classes(part->msh, g, class, 10, color);
}
