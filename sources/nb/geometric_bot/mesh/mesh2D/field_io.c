#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#include "nb/memory_bot.h"
#include "nb/io_bot.h"
#include "nb/geometric_bot.h"

#include "mesh2D_struct.h"

#define POW2(a) ((a)*(a)*(a)*(a))

#define STEP_READ_FINISH 100
#define FIELD_READ_FINISH 101
#define FIELD_SKIPPED 102

static int check_header(nb_cfreader_t *cfr);
static int check_mesh_correspondence(nb_cfreader_t *cfr,
				     const nb_mesh2D_t *mesh);

static int read_and_draw_all_steps(nb_cfreader_t *cfr,
				   const nb_mesh2D_t *mesh,
				   const char *dir_output,
				   void (*show_progress)(float prog));
static int read_step_data(nb_cfreader_t *cfr, const nb_mesh2D_t *mesh,
			  double *disp, double *damage,
			  uint32_t *step, double *time);
static int check_step_header(nb_cfreader_t *cfr, uint32_t *step,  double *time);
static int read_disp_field(nb_cfreader_t *cfr,
			   const nb_mesh2D_t *mesh,
			   double *disp);
static int read_vector(nb_cfreader_t *cfr, uint32_t id, double *vec);
static int read_damage_field(nb_cfreader_t *cfr,
			     const nb_mesh2D_t *mesh,
			     double *damage);
static int read_scalar(nb_cfreader_t *cfr, uint32_t id, double *scalar);
static void draw_damage(const char *dir_output, const nb_mesh2D_t *mesh,
			const double *disp, const double *damage,
			int step);
static void draw_nodal_damage_field(nb_graphics_context_t *g,
				    int width, int height,
				    const void *draw_data);
static void draw_face_damage_field(nb_graphics_context_t *g,
				    int width, int height,
				    const void *draw_data);
static int jump_mesh_correspondence(nb_cfreader_t *cfr, uint32_t *N_nodes,
				    uint32_t *N_edges, uint32_t *N_elems);
static int count_all_steps(nb_cfreader_t *cfr,
			   uint32_t N_nodes, uint32_t N_edges,
			   uint32_t N_elems, uint32_t *N_steps);
static int check_step(nb_cfreader_t *cfr,
		      uint32_t N_nodes, uint32_t N_edges,
		      uint32_t N_elems, uint32_t *readed_step);
static int check_field(nb_cfreader_t *cfr, uint32_t N_nodes,
		       uint32_t N_edges, uint32_t N_elems);
static int check_field_header(nb_cfreader_t *cfr, uint32_t N_nodes,
			      uint32_t N_edges, uint32_t N_elems,
			      int *N_comp, uint32_t *N, char *field_name);
static int get_N_components(const char *type);
static int get_field_length(const char *var, uint32_t N_nodes,
			    uint32_t N_edges, uint32_t N_elems);
static int jump_field_data(nb_cfreader_t *cfr, uint32_t N);
static int jump_to_step(nb_cfreader_t *cfr, uint32_t step_id,
			uint32_t N_nodes, uint32_t N_edges, uint32_t N_elems);
static int read_step_field(nb_cfreader_t *cfr, uint32_t N_nodes,
			   uint32_t N_edges, uint32_t N_elems, double *time,
			   const char *field_name, double *field);
static int read_field(nb_cfreader_t *cfr, uint32_t N_nodes,
		      uint32_t N_edges, uint32_t N_elems,
		      const char *fieldname, double *field);
static int read_field_data(nb_cfreader_t *cfr, uint32_t N,
			   int N_comp, double *field);
static int read_symtensor(nb_cfreader_t *cfr, uint32_t id, double *symtensor);
static int read_tensor(nb_cfreader_t *cfr, uint32_t id, double *tensor);

int nb_mesh2D_field_read_and_draw(const char *dir_saved_results,
				  const char *dir_output,
				  /* Show progress can be NULL */
				  void (*show_progress)(float prog))
{
	char name[100];
	sprintf(name, "%s/mesh.nbt", dir_saved_results);
	nb_mesh2D_type mesh_type;
	int status = nb_mesh2D_read_type_nbt(name, &mesh_type);
	if (0 != status)
		goto EXIT;

	uint32_t mesh_memsize = nb_mesh2D_get_memsize(mesh_type);
	nb_mesh2D_t *mesh = nb_soft_allocate_mem(mesh_memsize);
	nb_mesh2D_init(mesh, mesh_type);
	status = nb_mesh2D_read_nbt(mesh, name);
	if (0 != status)
		goto FINISH_MESH;

	nb_cfreader_t *cfr = nb_cfreader_create();
	nb_cfreader_load_nbt_format(cfr);
	sprintf(name, "%s/results.nbt", dir_saved_results);
	status = nb_cfreader_open_file(cfr, name);
	if (0 != status)
		goto DESTROY_READER;

	status = check_header(cfr);
	if (0 != status)
		goto CLOSE_FILE;

	status = check_mesh_correspondence(cfr, mesh);
	if (0 != status)
		goto CLOSE_FILE;

	status = read_and_draw_all_steps(cfr, mesh, dir_output,
					 show_progress);
	if (0 != status)
		goto CLOSE_FILE;

CLOSE_FILE:
	nb_cfreader_close_file(cfr);
DESTROY_READER:
	nb_cfreader_destroy(cfr);
FINISH_MESH:
	nb_mesh2D_finish(mesh);
	nb_soft_free_mem(mesh_memsize, mesh);
EXIT:
	return status;
}

static int check_header(nb_cfreader_t *cfr)
{
	char class[50];
	int status = nb_cfreader_nbt_check_header(cfr, class);
	if (0 != status)
		goto EXIT;

	if (0 != strcmp(class, "mesh2D_field"))
		status = 1;
EXIT:
	return status;
}

static int check_mesh_correspondence(nb_cfreader_t *cfr,
				     const nb_mesh2D_t *mesh)
{
	char var[50];
	int status = nb_cfreader_read_var_token(cfr, "Type", var);
	if (0 != status)
		goto EXIT;

	const char *mesh_type = nb_mesh2D_get_type_string(mesh);
	if (0 != strcmp(var, mesh_type)) {
		status = 1;
		goto EXIT;
	}
	
	uint32_t N;
	status = nb_cfreader_read_var_uint(cfr, "N_nodes", &N);
	if (0 != status)
		goto EXIT;

	if (N != nb_mesh2D_get_N_nodes(mesh)) {
		status = 1;
		goto EXIT;
	}

	status = nb_cfreader_read_var_uint(cfr, "N_edges", &N);
	if (0 != status)
		goto EXIT;

	if (N != nb_mesh2D_get_N_edges(mesh)) {
		status = 1;
		goto EXIT;
	}

	status = nb_cfreader_read_var_uint(cfr, "N_elems", &N);
	if (0 != status)
		goto EXIT;

	if (N != nb_mesh2D_get_N_elems(mesh)) {
		status = 1;
		goto EXIT;
	}
	status = 0;
EXIT:
	return status;

}

static int read_and_draw_all_steps(nb_cfreader_t *cfr,
				   const nb_mesh2D_t *mesh,
				   const char *dir_output,
				   void (*show_progress)(float prog))
{
	int status;
	uint32_t N = nb_mesh2D_get_N_elems(mesh);
	uint32_t memsize = 3 * N * sizeof(double);
	char *memblock = nb_allocate_mem(memsize);
	double *disp = (void*) memblock;
	double *damage = (void*) (memblock + 2 * N * sizeof(*disp));
	uint32_t step = 1;
	while (1) {
		uint32_t readed_step;
		double time;
		int read_status = read_step_data(cfr, mesh, disp, damage,
						 &readed_step, &time);
		if (STEP_READ_FINISH == read_status) {
			status = 0;
			goto EXIT;
		} else if (0 != read_status) {
			status = 1;
			goto EXIT;
		}
		if (readed_step != step) {
			status = 1;
			goto EXIT;
		}

		draw_damage(dir_output, mesh, disp, damage, step);
		step += 1;
		if (NULL != show_progress)
			show_progress(time);
	}
EXIT:
	nb_free_mem(memblock);
	return status;
}

static int read_step_data(nb_cfreader_t *cfr, const nb_mesh2D_t *mesh,
			  double *disp, double *damage, uint32_t *step,
			  double *time)
{
	int status = check_step_header(cfr, step, time);
	if (0 != status)
		goto EXIT;

	status = read_disp_field(cfr, mesh, disp);
	if (0 != status)
		goto EXIT;

	status = read_damage_field(cfr, mesh, damage);
	if (0 != status)
		goto EXIT;

	status = nb_cfreader_check_line(cfr, "End Step");
	if (0 != status)
		goto EXIT;
EXIT:
	return status;
}

static int check_step_header(nb_cfreader_t *cfr, uint32_t *step, double *time)
{
	int status = nb_cfreader_check_token(cfr, "Step");
	if (0 != status) {
		if (NB_CFREADER_EOF == status)
			status = STEP_READ_FINISH;
		else
			status = 1;
		goto EXIT;
	}
	
	status = nb_cfreader_read_uint(cfr, step);
	if (0 != status)
		goto EXIT;
	
	status = nb_cfreader_read_var_double(cfr, "Time", time);
	if (0 != status)
		goto EXIT;
EXIT:
	return status;
}

static int read_disp_field(nb_cfreader_t *cfr,
			   const nb_mesh2D_t *mesh,
			   double *disp)
{
	int N_comp;
	uint32_t N;
	char field_name[30];
	int status = check_field_header(cfr,
					nb_mesh2D_get_N_nodes(mesh),
					nb_mesh2D_get_N_edges(mesh),
					nb_mesh2D_get_N_elems(mesh),
					&N_comp, &N, field_name);
	if (0 != status)
		goto EXIT;

	if (N_comp != 2) {
		status = 1;
		goto EXIT;
	}

	if (nb_mesh2D_get_N_elems(mesh) != N) {
		status = 1;
		goto EXIT;
	}

	if (0 != strcmp("Displacement", field_name)) {
		status = 1;
		goto EXIT;
	}

	status = nb_cfreader_check_line(cfr, "Data");
	if (0 != status)
		goto EXIT;

	for (uint32_t i = 0; i < N; i++) {
		status = read_vector(cfr, i, &(disp[i*2]));
		if (0 != status)
			goto EXIT;
	}

	status = nb_cfreader_check_line(cfr, "End Data");
	if (0 != status)
		goto EXIT;

	status = nb_cfreader_check_line(cfr, "End Field");
	if (0 != status)
		goto EXIT;
EXIT:
	return status;
}

static int read_vector(nb_cfreader_t *cfr, uint32_t id, double *vec)
{
	uint32_t id_read;
	int status = nb_cfreader_read_uint(cfr, &id_read);
	if (0 != status)
		goto EXIT;
	if (id_read != id + 1) {
		status = 1;
		goto EXIT;
	}
	
	status = nb_cfreader_check_token(cfr, "<-");
	if (0 != status)
		goto EXIT;

	status = nb_cfreader_check_token(cfr, "[");
	if (0 != status)
		goto EXIT;

	double val;
	status = nb_cfreader_read_double(cfr, &val);
	if (0 != status)
		goto EXIT;
	vec[0] = val;

	status = nb_cfreader_check_token(cfr, ",");
	if (0 != status)
		goto EXIT;

	status = nb_cfreader_read_double(cfr, &val);
	if (0 != status)
		goto EXIT;
	vec[1] = val;

	status = nb_cfreader_check_token(cfr, "]");
	if (0 != status)
		goto EXIT;
EXIT:
	return status;
}

static int read_damage_field(nb_cfreader_t *cfr,
			     const nb_mesh2D_t *mesh,
			     double *damage)
{
	int N_comp;
	uint32_t N;
	char field_name[30];
	int status = check_field_header(cfr,
					nb_mesh2D_get_N_nodes(mesh),
					nb_mesh2D_get_N_edges(mesh),
					nb_mesh2D_get_N_elems(mesh),
					&N_comp, &N, field_name);
	if (0 != status)
		goto EXIT;

	if (N_comp != 1) {
		status = 1;
		goto EXIT;
	}

	if (nb_mesh2D_get_N_elems(mesh) != N) {
		status = 1;
		goto EXIT;
	}

	if (0 != strcmp("Damage", field_name)) {
		status = 1;
		goto EXIT;
	}

	status = nb_cfreader_check_line(cfr, "Data");
	if (0 != status)
		goto EXIT;
	
	for (uint32_t i = 0; i < N; i++) {
		status = read_scalar(cfr, i, &(damage[i]));
		if (0 != status)
			goto EXIT;
	}

	status = nb_cfreader_check_line(cfr, "End Data");
	if (0 != status)
		goto EXIT;

	status = nb_cfreader_check_line(cfr, "End Field");
	if (0 != status)
		goto EXIT;
EXIT:
	return status;
}

static int read_scalar(nb_cfreader_t *cfr, uint32_t id, double *scalar)
{
	uint32_t id_read;
	int status = nb_cfreader_read_uint(cfr, &id_read);
	if (0 != status)
		goto EXIT;
	if (id_read != id + 1) {
		status = 1;
		goto EXIT;
	}
	
	status = nb_cfreader_check_token(cfr, "<-");
	if (0 != status)
		goto EXIT;
	double val;
	status = nb_cfreader_read_double(cfr, &val);
	if (0 != status)
		goto EXIT;

	*scalar = val;

EXIT:
	return status;
}

static void draw_damage(const char *dir_output, const nb_mesh2D_t *mesh,
			const double *disp, const double *damage, int step)
{
	uint32_t N_nodes = nb_mesh2D_get_N_nodes(mesh);
	double *node_disp = nb_allocate_mem(2 * N_nodes * sizeof(*node_disp));
	double *node_damage = nb_allocate_mem(N_nodes * sizeof(*node_damage));
	nb_mesh2D_extrapolate_elems_to_nodes(mesh, 1, damage, node_damage);
	nb_mesh2D_extrapolate_elems_to_nodes(mesh, 2, disp, node_disp);
	
	char name[100];
	sprintf(name, "%s/damage_field_%i.png", dir_output, step);

	const void *data[4];
	data[0] = mesh;
	data[1] = node_damage;
	data[2] = disp;	
	data[3] = node_disp;	
	nb_graphics_export(name, 1000, 700, draw_nodal_damage_field, data);

	nb_free_mem(node_disp);
	nb_free_mem(node_damage);
}

static void get_centroid_color(const nb_mesh2D_t *msh,
			       const nb_palette_t *pal,
			       const double *field, uint32_t elem_id,
			       uint8_t cc[4])
/* Stolen from mshpoly_draw.c */
{
	uint16_t rgba_sum[4] = {0, 0, 0, 0};
	uint16_t N_adj = nb_mesh2D_elem_get_N_adj(msh, elem_id);
	for (uint16_t j = 0; j < N_adj; j++) {
		uint32_t nid = nb_mesh2D_elem_get_adj(msh, elem_id, j);
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

static void crack_distortion_field_on_nodes(const nb_mesh2D_t *msh,
					    const double *disp,
					    const double *nodal_disp,
					    const float disp_factor,
					    nb_graphics_context_t *g,
					    const double *normalized_field,
					    nb_palette_preset palette)
{
	nb_palette_t *pal = 
		nb_palette_create_preset(palette);

	uint32_t N_elems = nb_mesh2D_get_N_elems(msh);
	for (uint32_t i = 0; i < N_elems; i++) {
		double xc = nb_mesh2D_elem_get_x(msh, i);
		double yc = nb_mesh2D_elem_get_y(msh, i);
		double uxc = disp_factor * disp[i * 2];
		double uyc = disp_factor * disp[i*2+1];
		uint8_t cc[4];
		get_centroid_color(msh, pal, normalized_field, i, cc);

		uint16_t N_adj = nb_mesh2D_elem_get_N_adj(msh, i);
		for (uint16_t j = 0; j < N_adj; j++) {
			uint32_t n1 = nb_mesh2D_elem_get_adj(msh, i, j);
			double x1 = nb_mesh2D_node_get_x(msh, n1);
			double y1 = nb_mesh2D_node_get_y(msh, n1);
			double ux1 = disp_factor * nodal_disp[n1 * 2];
			double uy1 = disp_factor * nodal_disp[n1*2+1];
			double udist = sqrt(POW2((x1+ux1)-(xc+uxc)) +
					    POW2((y1+uy1)-(yc+uyc)));
			double xdist = sqrt(POW2(x1-xc) + POW2(y1-yc));
			if (udist/xdist > 2) {
				x1 += uxc;
				y1 += uyc;
			} else {
				x1 += ux1;
				y1 += uy1;
			}

			uint32_t n2 = nb_mesh2D_elem_get_adj(msh, i,
							      (j+1) % N_adj);
			double x2 = nb_mesh2D_node_get_x(msh, n2);
			double y2 = nb_mesh2D_node_get_y(msh, n2);
			double ux2 = disp_factor * nodal_disp[n2 * 2];
			double uy2 = disp_factor * nodal_disp[n2*2+1];
			udist = sqrt(POW2((x2+ux2)-(xc+uxc)) +
					    POW2((y2+uy2)-(yc+uyc)));
			xdist = sqrt(POW2(x2-xc) + POW2(y2-yc));
			if (udist/xdist > 2) {
				x2 += uxc;
				y2 += uyc;
			} else {
				x2 += ux2;
				y2 += uy2;
			}

			nb_graphics_move_to(g, xc + uxc, yc + uyc);
			nb_graphics_line_to(g, x1, y1);
			nb_graphics_line_to(g, x2, y2);
			nb_graphics_close_path(g);
			
			uint8_t c1[4], c2[4];
			nb_palette_get_rgba(pal, normalized_field[n1], c1);
			nb_palette_get_rgba(pal, normalized_field[n2], c2);

			nb_graphics_set_source_trg(g, xc + uxc, yc + uyc,
						   x1, y1, x2, y2,
						   cc, c1, c2);
			nb_graphics_fill(g);
		}
	}
	nb_palette_destroy(pal);
}

static void draw_nodal_damage_field(nb_graphics_context_t *g,
				    int width, int height,
				    const void *draw_data)
{
	void **data = (void*) draw_data;
	const nb_mesh2D_t *mesh = data[0];
	const double *damage = data[1];
	const double *disp = data[2];
	const double *nodal_disp = data[3];

	if (!nb_graphics_is_camera_enabled(g)) {
		double box[4];
		mesh->get_enveloping_box(mesh->msh, box);

		nb_graphics_enable_camera(g);
		nb_graphics_camera_t* cam = nb_graphics_get_camera(g);
		nb_graphics_cam_fit_box(cam, box, width, height);
	}

	nb_graphics_disable_camera(g);
	nb_graphics_set_source(g, NB_WHITE);
	nb_graphics_set_rectangle(g, 0, 0, width, height);
	nb_graphics_fill(g);
	nb_graphics_enable_camera(g);
	
	crack_distortion_field_on_nodes(mesh, disp, nodal_disp, 1e1,
					g, damage, NB_RAINBOW);

}

static void draw_face_damage_field(nb_graphics_context_t *g,
				    int width, int height,
				    const void *draw_data)
{
	void **data = (void*) draw_data;
	const nb_mesh2D_t *mesh = data[0];
	const double *damage = data[1];

	if (!nb_graphics_is_camera_enabled(g)) {
		double box[4];
		mesh->get_enveloping_box(mesh->msh, box);

		nb_graphics_enable_camera(g);
		nb_graphics_camera_t* cam = nb_graphics_get_camera(g);
		nb_graphics_cam_fit_box(cam, box, width, height);
	}

	nb_graphics_disable_camera(g);
	nb_graphics_set_source(g, NB_WHITE);
	nb_graphics_set_rectangle(g, 0, 0, width, height);
	nb_graphics_fill(g);
	nb_graphics_enable_camera(g);

	nb_graphics_set_source(g, NB_BLACK);
	mesh->graphics.fill_elems(mesh->msh, g);
	mesh->graphics.draw_field_on_faces(mesh->msh, g,
					   damage, NB_SUNSET);

	nb_graphics_set_source(g, NB_BLACK);
	nb_graphics_set_line_width(g, 1.5);
	mesh->graphics.draw_boundaries(mesh->msh, g);
}

int nb_mesh2D_field_get_N_steps(const char *filename, uint32_t *N_steps)
{
	nb_cfreader_t *cfr = nb_cfreader_create();
	nb_cfreader_load_nbt_format(cfr);
	
	int status = nb_cfreader_open_file(cfr, filename);
	if (0 != status)
		goto EXIT;

	status = check_header(cfr);
	if (0 != status)
		goto CLOSE_FILE;

	uint32_t N_nodes;
	uint32_t N_edges;
	uint32_t N_elems;
	status = jump_mesh_correspondence(cfr, &N_nodes, &N_edges, &N_elems);
	if (0 != status)
		goto CLOSE_FILE;

	uint32_t N;
	status = count_all_steps(cfr, N_nodes, N_edges, N_elems, &N);
	if (0 != status)
		goto CLOSE_FILE;

	status = 0;
	*N_steps = N;
CLOSE_FILE:
	nb_cfreader_close_file(cfr);
EXIT:
	nb_cfreader_destroy(cfr);
	return status;
}

static int jump_mesh_correspondence(nb_cfreader_t *cfr, uint32_t *N_nodes,
				    uint32_t *N_edges, uint32_t *N_elems)
{
	char var[50];
	int status = nb_cfreader_read_var_token(cfr, "Type", var);
	if (0 != status)
		goto EXIT;
	
	uint32_t N;
	status = nb_cfreader_read_var_uint(cfr, "N_nodes", &N);
	if (0 != status)
		goto EXIT;
	*N_nodes = N;

	status = nb_cfreader_read_var_uint(cfr, "N_edges", &N);
	if (0 != status)
		goto EXIT;
	*N_edges = N;

	status = nb_cfreader_read_var_uint(cfr, "N_elems", &N);
	if (0 != status)
		goto EXIT;
	*N_elems = N;
	status = 0;
EXIT:
	return status;

}

static int count_all_steps(nb_cfreader_t *cfr,
			   uint32_t N_nodes, uint32_t N_edges,
			   uint32_t N_elems, uint32_t *N_steps)
{
	int status;
	int step = 1;
	while (1) {
		uint32_t readed_step;
		int read_status = check_step(cfr, N_nodes, N_edges,
					     N_elems, &readed_step);
		if (STEP_READ_FINISH == read_status) {
			status = 0;
			*N_steps = step - 1;
			goto EXIT;
		} else if (0 != read_status) {
			status = 1;
			goto EXIT;
		}
		if (readed_step != step) {
			status = 1;
			goto EXIT;
		}
		step += 1;
	}
EXIT:
	return status;
}

static int check_step(nb_cfreader_t *cfr,
		      uint32_t N_nodes, uint32_t N_edges,
		      uint32_t N_elems, uint32_t *readed_step)
{
	double time;
	uint32_t step;
	int status = check_step_header(cfr, &step, &time);
	if (0 != status)
		goto EXIT;

	while (1) {
		status = check_field(cfr, N_nodes, N_edges, N_elems);
		if (FIELD_READ_FINISH == status)
			break;
		if (0 != status)
			goto EXIT;
	}

	status = nb_cfreader_check_line(cfr, "End Step");
	if (0 != status)
		goto EXIT;

	*readed_step = step;
EXIT:
	return status;
}

static int check_field(nb_cfreader_t *cfr, uint32_t N_nodes,
		       uint32_t N_edges, uint32_t N_elems)
{
	int N_comp;
	uint32_t N;
	char field_name[30];
	int status = check_field_header(cfr, N_nodes, N_edges, N_elems,
					&N_comp, &N, field_name);
	if (0 != status)
		goto EXIT;

	status = nb_cfreader_check_line(cfr, "Data");
	if (0 != status)
		goto EXIT;

	status = jump_field_data(cfr, N);
	if (0 != status)
		goto EXIT;

	status = nb_cfreader_check_line(cfr, "End Data");
	if (0 != status)
		goto EXIT;

	status = nb_cfreader_check_line(cfr, "End Field");
	if (0 != status)
		goto EXIT;
EXIT:
	return status;
}

static int check_field_header(nb_cfreader_t *cfr, uint32_t N_nodes,
			      uint32_t N_edges, uint32_t N_elems,
			      int *N_comp, uint32_t *N, char *field_name)
{
	int status = nb_cfreader_check_line(cfr, "Field");
	if (0 != status) {
		status = FIELD_READ_FINISH;
		goto EXIT;
	}
	
	char var[30];
	status = nb_cfreader_read_var_token(cfr, "Type", var);
	if (0 != status)
		goto EXIT;

	*N_comp = get_N_components(var);

	if (0 == *N_comp) {
		status = 1;
		goto EXIT;
	}

	status = nb_cfreader_read_var_string(cfr, "Name", field_name);
	if (0 != status)
		goto EXIT;

	status = nb_cfreader_read_var_token(cfr, "Support", var);
	if (0 != status)
		goto EXIT;

	*N = get_field_length(var, N_nodes, N_edges, N_elems);
	if (0 == *N) {
		status = 1;
		goto EXIT;
	}
EXIT:
	return status;
}

static int get_N_components(const char *type)
{
	int out;
	if (0 == strcmp(type, "Scalar"))
		out = 1;
	else if (0 == strcmp(type, "Vector"))
		out = 2;
	else if (0 == strcmp(type, "SymTensor"))
		out = 3;
	else if (0 == strcmp(type, "Tensor"))
		out = 4;
	else
		out = 0;
	return out;
}

static int get_field_length(const char *var, uint32_t N_nodes,
			    uint32_t N_edges, uint32_t N_elems)
{
	int out;
	if (0 == strcmp(var, "Nodes"))
		out = N_nodes;
	else if (0 == strcmp(var, "Edges"))
		out = N_edges;
	else if (0 == strcmp(var, "Elements"))
		out = N_elems;
	else
		out = 0;
	return out;

}

static int jump_field_data(nb_cfreader_t *cfr, uint32_t N)
{
	int status = 0;
	for (uint32_t i = 0; i < N; i++) {
		status = nb_cfreader_jump_line(cfr);
		if (0 != status)
			goto EXIT;
	}
EXIT:
	return status;
}

int nb_mesh2D_field_read_step_field(const char *filename, uint32_t step,
				    double *time, const char *field_name,
				    double *field)
{
	nb_cfreader_t *cfr = nb_cfreader_create();
	nb_cfreader_load_nbt_format(cfr);
	
	int status = nb_cfreader_open_file(cfr, filename);
	if (0 != status)
		goto EXIT;

	status = check_header(cfr);
	if (0 != status)
		goto CLOSE_FILE;

	uint32_t N_nodes;
	uint32_t N_edges;
	uint32_t N_elems;
	status = jump_mesh_correspondence(cfr, &N_nodes, &N_edges, &N_elems);
	if (0 != status)
		goto CLOSE_FILE;

	status = jump_to_step(cfr, step, N_nodes, N_edges, N_elems);
	if (0 != status)
		goto CLOSE_FILE;

	status = read_step_field(cfr, N_nodes, N_edges, N_elems, time,
				 field_name, field);
	if (0 != status)
		goto CLOSE_FILE;

	status = 0;
CLOSE_FILE:
	nb_cfreader_close_file(cfr);
EXIT:
	nb_cfreader_destroy(cfr);
	return status;
}

static int jump_to_step(nb_cfreader_t *cfr, uint32_t step_id,
			uint32_t N_nodes, uint32_t N_edges, uint32_t N_elems)
{
	int status;
	for (uint32_t i = 1; i < step_id; i++) {
		uint32_t readed_step;
		int read_status = check_step(cfr, N_nodes, N_edges,
					     N_elems, &readed_step);
		if (0 != read_status) {
			status = 1;
			goto EXIT;
		}
		if (readed_step != i) {
			status = 1;
			goto EXIT;
		}
	}
EXIT:
	return status;
}

static int read_step_field(nb_cfreader_t *cfr, uint32_t N_nodes,
			   uint32_t N_edges, uint32_t N_elems, double *time,
			   const char *field_name, double *field)
{
	uint32_t step;
	int status = check_step_header(cfr, &step, time);
	if (0 != status)
		goto EXIT;

	bool found = false;
	while (1) {
		status = read_field(cfr, N_nodes, N_edges, N_elems,
				    field_name, field);
		if (0 == status)
			found = true;
		if (FIELD_READ_FINISH == status)
			break;
		if (0 != status && FIELD_SKIPPED != status) {
			status = 1;
			goto EXIT;
		}
	}

	status = nb_cfreader_check_line(cfr, "End Step");
	if (0 != status)
		goto EXIT;

	if (!found)
		status = 1;
EXIT:
	return status;	
}

static int read_field(nb_cfreader_t *cfr, uint32_t N_nodes,
		      uint32_t N_edges, uint32_t N_elems,
		      const char *fieldname, double *field)
{
	int N_comp;
	uint32_t N;
	char fieldname_readed[30];
	int status = check_field_header(cfr, N_nodes, N_edges, N_elems,
					&N_comp, &N, fieldname_readed);
	if (0 != status)
		goto EXIT;

	status = nb_cfreader_check_line(cfr, "Data");
	if (0 != status)
		goto EXIT;
	
	int field_match = strcmp(fieldname_readed, fieldname);
	if (0 != field_match) {
		status = jump_field_data(cfr, N);
		if (0 != status)
			goto EXIT;
	} else {
		status = read_field_data(cfr, N, N_comp, field);
		if (0 != status)
			goto EXIT;
	}

	status = nb_cfreader_check_line(cfr, "End Data");
	if (0 != status)
		goto EXIT;

	status = nb_cfreader_check_line(cfr, "End Field");
	if (0 != status)
		goto EXIT;
	if (0 != field_match)
		status = FIELD_SKIPPED;
EXIT:
	return status;
}

static int read_field_data(nb_cfreader_t *cfr, uint32_t N,
			   int N_comp, double *field)
{
	int status;
	int (*read_val)(nb_cfreader_t *cfr, uint32_t id, double *val);
	if (1 == N_comp) {
		read_val = read_scalar;
	} else if (2 == N_comp) {
		read_val = read_vector;
	} else if (3 == N_comp) {
		read_val = read_symtensor;
	} else if (4 == N_comp) {
		read_val = read_tensor;
	} else {
		status = 1;
		goto EXIT;
	}
	
	for (uint32_t i = 0; i < N; i++) {
		status = read_val(cfr, i, &(field[i * N_comp]));
		if (0 != status)
			goto EXIT;
	}
EXIT:
	return status;
}

static int read_symtensor(nb_cfreader_t *cfr, uint32_t id, double *symtensor)
{
	return 1;/* TEMPORAL */
}

static int read_tensor(nb_cfreader_t *cfr, uint32_t id, double *tensor)
{
	return 1;/* TEMPORAL */
}
