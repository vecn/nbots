#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>

#include "nb/memory_bot.h"
#include "nb/container_bot.h"
#include "nb/graphics_bot.h"
#include "nb/geometric_bot.h"

#include "mesh2D_struct.h"

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
#define COLOR_LEVEL_SET NB_BLACK

#define PALETTE_W 30
#define PALETTE_H 250
#define PALETTE_MARGIN 40

typedef struct {
	const nb_mesh2D_t *mesh;
	const void *values;
	nb_mesh2D_entity vals_entity;
	nb_mesh2D_array_type vals_type;
	bool draw_wires;
} draw_data;

typedef struct {
	const nb_mesh2D_t *mesh;
	const double *field;
	uint16_t N_ls;
	bool draw_wires;
} level_set_data;

static void init_draw_data(draw_data *data,
			   const nb_mesh2D_t *mesh,
			   nb_mesh2D_entity vals_entity,
			   nb_mesh2D_array_type vals_type,
			   const void *values,
			   bool draw_wires);

static void draw(nb_graphics_context_t *g, int width, int height,
		 const void *draw_data);
static void set_camera(nb_graphics_context_t *g, int width, int height,
		       const nb_mesh2D_t *mesh);
static void fill(const nb_mesh2D_t *mesh,
		 nb_graphics_context_t *g,
		 const draw_data *data);

static void fill_elems_field_on_nodes(const nb_mesh2D_t *mesh,
				      nb_graphics_context_t *g,
				      const double *values);
static void normalize_values(double *normalized_values,
			     const double *values, uint32_t N);
static void fill_elems_field_on_elems(const nb_mesh2D_t *mesh,
				      nb_graphics_context_t *g,
				      const double *values);
static void fill_elems_classes(const nb_mesh2D_t *mesh,
			       nb_graphics_context_t *g,
			       const uint8_t *class);
static void set_class_colors(nb_graphics_color_t color[10]);
static void fill_elems(const nb_mesh2D_t *mesh,
		       nb_graphics_context_t *g);

static void draw_faces(const nb_mesh2D_t *mesh,
		       nb_graphics_context_t *g,
		       const draw_data *data);

static void draw_field_on_faces(const nb_mesh2D_t *mesh,
				nb_graphics_context_t *g,
				const double *values);
static void draw_classes_on_faces(const nb_mesh2D_t *mesh,
				  nb_graphics_context_t *g,
				  const uint8_t *class);
static void draw_wires(const nb_mesh2D_t *mesh,
		       nb_graphics_context_t *g);

static void draw_boundaries(const nb_mesh2D_t *mesh,
			    nb_graphics_context_t *g);

static void fill_nodes_classes(const nb_mesh2D_t *mesh,
			       nb_graphics_context_t *g,
			       const uint8_t *class);
static void add_palette(const nb_mesh2D_t *mesh,
			const draw_data *data,
			nb_graphics_context_t *g, int width, int height);
static void get_min_max(const nb_mesh2D_t *mesh,
			const draw_data *data,
			double *min, double *max);
static void init_level_set_data(level_set_data *data,
				const nb_mesh2D_t *mesh,
				const double *field, uint16_t N_ls,
				bool draw_wires);
static void draw_level_sets(nb_graphics_context_t *g, int width, int height,
			    const void *draw_data);
static void process_level_sets(const level_set_data *data,
			       nb_graphics_context_t *g);

void nb_mesh2D_export_draw(const nb_mesh2D_t *mesh,
			      const char *filename,
			      int width, int height,
			      nb_mesh2D_entity vals_entity,
			      nb_mesh2D_array_type vals_type,
			      const void *values,
			      bool draw_wires)
{
	draw_data data;
	init_draw_data(&data, mesh, vals_entity,
		       vals_type, values, draw_wires);

	nb_graphics_export(filename, width, height, draw, &data);
}

static void init_draw_data(draw_data *data,
			   const nb_mesh2D_t *mesh,
			   nb_mesh2D_entity vals_entity,
			   nb_mesh2D_array_type vals_type,
			   const void *values,
			   bool draw_wires)
{
	data->mesh = mesh;
	data->values = values;
	data->vals_entity = vals_entity;
	data->vals_type = vals_type;
	data->draw_wires = draw_wires;
}

static void draw(nb_graphics_context_t *g, int width, int height,
		 const void *data_ptr)
{
	const draw_data *data = data_ptr;
	const nb_mesh2D_t *mesh = data->mesh;

	if (!nb_graphics_is_camera_enabled(g))
		set_camera(g, width, height, mesh);

	fill(mesh, g, data);

	draw_faces(mesh, g, data);

	draw_boundaries(mesh, g);

	if (NB_NODE == data->vals_entity && NB_CLASS == data->vals_type)
		fill_nodes_classes(mesh, g, data->values);

	if (NB_FIELD == data->vals_type)
		add_palette(mesh, data, g, width, height);
}

static void set_camera(nb_graphics_context_t *g, int width, int height,
		       const nb_mesh2D_t *mesh)
{
	double box[4];
	mesh->get_enveloping_box(mesh->msh, box);

	nb_graphics_enable_camera(g);
	nb_graphics_camera_t* cam = nb_graphics_get_camera(g);
	nb_graphics_cam_fit_box(cam, box, width, height);
}

static void fill(const nb_mesh2D_t *mesh,
		 nb_graphics_context_t *g,
		 const draw_data *data)
{
	nb_mesh2D_entity enty = data->vals_entity;
	nb_mesh2D_array_type type = data->vals_type;

	if (NB_NODE == enty && NB_FIELD == type)
		fill_elems_field_on_nodes(mesh, g, data->values);

	else if (NB_ELEMENT == enty && NB_FIELD == type)
		fill_elems_field_on_elems(mesh, g, data->values);

	else if (NB_ELEMENT == enty && NB_CLASS == type)
		fill_elems_classes(mesh, g, data->values);
	else
		fill_elems(mesh, g);
}

static void fill_elems_field_on_nodes(const nb_mesh2D_t *mesh,
				      nb_graphics_context_t *g,
				      const double *values)
{
	uint32_t N = mesh->get_N_nodes(mesh->msh);
	uint32_t memsize = N * sizeof(double);
	double *normalized_values = nb_soft_allocate_mem(memsize);

	normalize_values(normalized_values, values, N);
	
	mesh->graphics.fill_elems_field_on_nodes(mesh->msh, g,
						 normalized_values,
						 PALETTE_FIELD);

	nb_soft_free_mem(memsize, normalized_values);
}

static void normalize_values(double *normalized_values,
			     const double *values, uint32_t N)
{
	uint32_t min_id;
	uint32_t max_id;
	nb_array_get_min_max_ids(values, N, sizeof(*values),
				  nb_compare_double,
				  &min_id, &max_id);

	double min = values[min_id];
	double range = values[max_id] - min;
	for (uint32_t i = 0; i < N; i++)
		normalized_values[i] = (values[i] - min) / range;
}

static void fill_elems_field_on_elems(const nb_mesh2D_t *mesh,
				      nb_graphics_context_t *g,
				      const double *values)
{
	uint32_t N = mesh->get_N_elems(mesh->msh);
	uint32_t memsize = N * sizeof(double);
	double *normalized_values = nb_soft_allocate_mem(memsize);

	normalize_values(normalized_values, values, N);
	
	mesh->graphics.fill_elems_field_on_elems(mesh->msh, g,
						 normalized_values,
						 PALETTE_FIELD);

	nb_soft_free_mem(memsize, normalized_values);
}

static void fill_elems_classes(const nb_mesh2D_t *mesh,
			       nb_graphics_context_t *g,
			       const uint8_t *class)
{
	nb_graphics_color_t color[10];
	set_class_colors(color);
	mesh->graphics.fill_elems_classes(mesh->msh, g, class, 10, color);
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

static void fill_elems(const nb_mesh2D_t *mesh,
		       nb_graphics_context_t *g)
{
	nb_graphics_set_source(g, COLOR_ELEM);
	mesh->graphics.fill_elems(mesh->msh, g);
}

static void draw_faces(const nb_mesh2D_t *mesh,
		       nb_graphics_context_t *g,
		       const draw_data *data)
{
	nb_mesh2D_entity enty = data->vals_entity;
	nb_mesh2D_array_type type = data->vals_type;
	if (NB_FACE == enty && NB_FIELD == type)
		draw_field_on_faces(mesh, g, data->values);
	else if (NB_FACE == enty && NB_CLASS == type)
		draw_classes_on_faces(mesh, g, data->values);
	else if (data->draw_wires)
		draw_wires(mesh, g);
}

static void draw_field_on_faces(const nb_mesh2D_t *mesh,
				nb_graphics_context_t *g,
				const double *values)
{
	nb_graphics_set_line_width(g, 1.5);
	mesh->graphics.draw_field_on_faces(mesh->msh, g, values,
					   PALETTE_FIELD);
}

static void draw_classes_on_faces(const nb_mesh2D_t *mesh,
				  nb_graphics_context_t *g,
				  const uint8_t *class)
{
	nb_graphics_color_t color[10];
	set_class_colors(color);
	nb_graphics_set_line_width(g, 1.5);
	mesh->graphics.draw_classes_on_faces(mesh->msh, g, class, 10, color);
}

static void draw_wires(const nb_mesh2D_t *mesh,
		       nb_graphics_context_t *g)
{
	nb_graphics_set_source(g, COLOR_EDGE);
	nb_graphics_set_line_width(g, 1.0);
	mesh->graphics.draw_wires(mesh->msh, g);
}

static void draw_boundaries(const nb_mesh2D_t *mesh,
			    nb_graphics_context_t *g)
{
	nb_graphics_set_source(g, COLOR_SGM);
	nb_graphics_set_line_width(g, 1.5);
	mesh->graphics.draw_boundaries(mesh->msh, g);
}

static void fill_nodes_classes(const nb_mesh2D_t *mesh,
			       nb_graphics_context_t *g,
			       const uint8_t *class)
{
	nb_graphics_color_t color[10];
	set_class_colors(color);
	mesh->graphics.fill_nodes_classes(mesh->msh, g, class, 10, color);
}

static void add_palette(const nb_mesh2D_t *mesh,
			const draw_data *data,
			nb_graphics_context_t *g, int width, int height)
{
	nb_palette_t *palette =
		nb_palette_create_preset(PALETTE_FIELD);
	double min, max;
	get_min_max(mesh, data, &min, &max);
	
	float label_width = 55;
	nb_palette_draw(g, palette,
				 width - PALETTE_W - PALETTE_MARGIN -
				 label_width,
				 height - PALETTE_H - PALETTE_MARGIN,
				 PALETTE_W,
				 PALETTE_H,
				 1, min, max);
	nb_palette_destroy(palette);
}

static void get_min_max(const nb_mesh2D_t *mesh,
			const draw_data *data,
			double *min, double *max)
{
	uint32_t N;
	if (NB_ELEMENT == data->vals_entity)
		N = nb_mesh2D_get_N_elems(mesh);
	else
		N = nb_mesh2D_get_N_nodes(mesh);

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

void nb_mesh2D_export_level_sets(const nb_mesh2D_t *mesh,
				 const char *filename, int width, int height,
				 const void *field, uint16_t N_ls,
				 bool draw_wires)
{
	level_set_data data;
	init_level_set_data(&data, mesh, field, N_ls, draw_wires);

	nb_graphics_export(filename, width, height, draw_level_sets, &data);
}

static void init_level_set_data(level_set_data *data,
				const nb_mesh2D_t *mesh,
				const double *field, uint16_t N_ls,
				bool draw_wires)
{
	data->mesh = mesh;
	data->field = field;
	data->N_ls = N_ls;
	data->draw_wires = draw_wires;
}

static void draw_level_sets(nb_graphics_context_t *g, int width, int height,
			    const void *draw_data)
{
	const level_set_data *data = draw_data;
	const nb_mesh2D_t *mesh = data->mesh;

	if (!nb_graphics_is_camera_enabled(g))
		set_camera(g, width, height, mesh);

	process_level_sets(data, g);

	if (data->draw_wires)
		draw_wires(mesh, g);

	draw_boundaries(mesh, g);
}

static void process_level_sets(const level_set_data *data,
			       nb_graphics_context_t *g)
{
	const nb_mesh2D_t *mesh = data->mesh;
	const double *field = data->field;

	nb_graphics_set_source(g, COLOR_LEVEL_SET);
	nb_graphics_set_line_width(g, 1.5);

	uint32_t min_id;
	uint32_t max_id;
	uint32_t N_nodes = nb_mesh2D_get_N_nodes(mesh);
	nb_array_get_min_max_ids(field, N_nodes, sizeof(double),
				 nb_compare_double, &min_id, &max_id);
	double min = field[min_id];
	double max = field[max_id];

	uint16_t N = data->N_ls;
	double level_set_step = (max - min) / (N - 1.0);
	for (uint16_t i = 0; i < N; i++) {
		double level_set = min + i * level_set_step;
		nb_mesh2D_draw_level_set(mesh, g, field, level_set);
	}
}

void nb_mesh2D_draw_wires(const nb_mesh2D_t *mesh,
			  nb_graphics_context_t *g)
{
	mesh->graphics.draw_wires(mesh->msh, g);
}

void nb_mesh2D_draw_boundaries(const nb_mesh2D_t *mesh,
				  nb_graphics_context_t *g)
{
	mesh->graphics.draw_boundaries(mesh->msh, g);
}

void nb_mesh2D_fill_elems(const nb_mesh2D_t *mesh,
			     nb_graphics_context_t *g)
{
	mesh->graphics.fill_elems(mesh->msh, g);
}

void nb_mesh2D_fill_elems_field_on_nodes(const nb_mesh2D_t *mesh,
					 nb_graphics_context_t *g,
					 const double *normalized_field,
					 nb_palette_preset palette)
{
	mesh->graphics.fill_elems_field_on_nodes(mesh->msh, g,
						 normalized_field,
						 palette);
}

void nb_mesh2D_fill_elems_field_on_elems(const nb_mesh2D_t *mesh,
					    nb_graphics_context_t *g,
					    const double *normalized_field,
					    nb_palette_preset palette)
{
	mesh->graphics.fill_elems_field_on_elems(mesh->msh, g,
						 normalized_field,
						 palette);
}

void nb_mesh2D_fill_elems_classes(const nb_mesh2D_t *mesh,
				     nb_graphics_context_t *g,
				     const uint8_t *class, uint8_t N_colors,
				     const nb_graphics_color_t *colors)
{
	mesh->graphics.fill_elems_classes(mesh->msh, g, class,
					  N_colors, colors);
}

void nb_mesh2D_fill_nodes(const nb_mesh2D_t *mesh,
			  nb_graphics_context_t *g)
{
	mesh->graphics.fill_nodes(mesh->msh, g);
}

void nb_mesh2D_fill_nodes_classes(const nb_mesh2D_t *mesh,
				  nb_graphics_context_t *g,
				  const uint8_t *class, uint8_t N_colors,
				  const nb_graphics_color_t *colors)
{
	mesh->graphics.fill_nodes_classes(mesh->msh, g, class,
					  N_colors, colors);
}

void nb_mesh2D_draw_field_on_faces(const nb_mesh2D_t *mesh,
				   nb_graphics_context_t *g,
				   const double *normalized_field,
				   nb_palette_preset palette)
{
	mesh->graphics.draw_field_on_faces(mesh->msh, g, normalized_field,
					   palette);
}

void nb_mesh2D_draw_classes_on_faces(const nb_mesh2D_t *mesh,
				     nb_graphics_context_t *g,
				     const uint8_t *class, uint8_t N_colors,
				     const nb_graphics_color_t *colors)
{
	mesh->graphics.draw_classes_on_faces(mesh->msh, g, class,
					     N_colors, colors);
}

void nb_mesh2D_draw_level_set(const nb_mesh2D_t *mesh,
			      nb_graphics_context_t *g,
			      const double *field_on_nodes,
			      double level_set)
{
	mesh->graphics.draw_level_set(mesh->msh, g, field_on_nodes,
				      level_set);
}
