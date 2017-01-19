#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>

#include "nb/memory_bot.h"
#include "nb/io_bot.h"
#include "nb/geometric_bot.h"

#define POW2(a) ((a)*(a)*(a)*(a))

#define STEP_READ_FINISH 100
#define FIELD_READ_FINISH 101

static int check_header(nb_cfreader_t *cfr);
static int check_mesh_correspondence(nb_cfreader_t *cfr,
				     const nb_mesh2D_t *mesh);

static int read_and_draw_all_steps(nb_cfreader_t *cfr,
				   const nb_mesh2D_t *mesh,
				   const char *dir_output,
				   void (*show_progress)(float prog));
static int read_step_data(nb_cfreader_t *cfr, const nb_mesh2D_t *mesh,
			  double *disp, double *damage,
			  int *step, double *time);
static int read_disp_field(nb_cfreader_t *cfr,
			   const nb_mesh2D_t *mesh,
			   double *disp);
static int read_vector(nb_cfreader_t *cfr, uint32_t id, double *vec);
static int read_damage_field(nb_cfreader_t *cfr,
			     const nb_mesh2D_t *mesh,
			     double *damage);
static int read_scalar(nb_cfreader_t *cfr, uint32_t id, double *scalar);
static void draw_damage(const char *dir_output, const nb_mesh2D_t *mesh,
			const double *damage, int step);
static int jump_mesh_correspondence(nb_cfreader_t *cfr, uint32_t *N_nodes,
				    uint32_t *N_edges, uint32_t *N_elems);
static int count_all_steps(nb_cfreader_t *cfr,
			   uint32_t N_nodes, uint32_t N_edges,
			   uint32_t N_elems, uint32_t *N_steps);
static int check_step(nb_cfreader_t *cfr,
		      uint32_t N_nodes, uint32_t N_edges,
		      uint32_t N_elems, int *readed_step);
static int check_field(nb_cfreader_t *cfr, uint32_t N_nodes,
		       uint32_t N_edges, uint32_t N_elems);
static int get_N_components(const char *type);
static int get_field_length(const char *var, uint32_t N_nodes,
			    uint32_t N_edges, uint32_t N_elems);

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
	int step = 1;
	while (1) {
		int readed_step;
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

		draw_damage(dir_output, mesh, damage, step);
		step += 1;
		if (NULL != show_progress)
			show_progress(time);
	}
EXIT:
	nb_free_mem(memblock);
	return status;
}

static int read_step_data(nb_cfreader_t *cfr, const nb_mesh2D_t *mesh,
			  double *disp, double *damage, int *step,
			  double *time)
{
	int status = nb_cfreader_check_token(cfr, "Step");
	if (0 != status) {
		if (NB_CFREADER_EOF == status)
			status = STEP_READ_FINISH;
		else
			status = 1;
		goto EXIT;
	}
	
	status = nb_cfreader_read_int(cfr, step);
	if (0 != status)
		goto EXIT;
	
	status = nb_cfreader_read_var_double(cfr, "Time", time);
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

static int read_disp_field(nb_cfreader_t *cfr,
			   const nb_mesh2D_t *mesh,
			   double *disp)
{
	int status = nb_cfreader_check_line(cfr, "Field");
	if (0 != status)
		goto EXIT;
	
	status = nb_cfreader_check_line(cfr, "Type = Vector");
	if (0 != status)
		goto EXIT;
	status = nb_cfreader_check_line(cfr, "Name = \"Displacement\"");
	if (0 != status)
		goto EXIT;
	status = nb_cfreader_check_line(cfr, "Support = Elements");
	if (0 != status)
		goto EXIT;
	status = nb_cfreader_check_line(cfr, "Data");
	if (0 != status)
		goto EXIT;

	uint32_t N = nb_mesh2D_get_N_elems(mesh);
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
	int status = nb_cfreader_check_line(cfr, "Field");
	if (0 != status)
		goto EXIT;
	
	status = nb_cfreader_check_line(cfr, "Type = Scalar");
	if (0 != status)
		goto EXIT;
	status = nb_cfreader_check_line(cfr, "Name = \"Damage\"");
	if (0 != status)
		goto EXIT;
	status = nb_cfreader_check_line(cfr, "Support = Elements");
	if (0 != status)
		goto EXIT;
	status = nb_cfreader_check_line(cfr, "Data");
	if (0 != status)
		goto EXIT;
	
	uint32_t N = nb_mesh2D_get_N_elems(mesh);
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
			const double *damage, int step)
{
	uint32_t N_nodes = nb_mesh2D_get_N_nodes(mesh);

	double *node_damage = nb_allocate_mem(N_nodes * sizeof(*node_damage));
	nb_mesh2D_extrapolate_elems_to_nodes(mesh, 1, damage, node_damage);
	
	char name[100];
	sprintf(name, "%s/damage_field_%i.png", dir_output, step);

	nb_mesh2D_export_draw(mesh, name, 1000, 700, NB_NODE, NB_FIELD,
			      node_damage, true);

	uint32_t N_faces = nb_mesh2D_get_N_edges(mesh);
	double *face_damage = nb_allocate_mem(N_faces * sizeof(double));

	for (uint32_t i = 0; i < N_faces; i++) {		
		uint32_t id1 = nb_mesh2D_edge_get_1n(mesh, i);
		uint32_t id2 = nb_mesh2D_edge_get_2n(mesh, i);
		double dmg1 = POW2(node_damage[id1]);
		double dmg2 = POW2(node_damage[id2]);
		face_damage[i] = (dmg1 + dmg2) / 2;
	}

	sprintf(name, "%s/fracture_%i.png", dir_output, step);

	nb_mesh2D_export_draw(mesh, name, 1000, 800, NB_FACE, NB_FIELD,
			      face_damage, true);

	nb_free_mem(face_damage);

	nb_free_mem(node_damage);
}

int nb_mesh2D_field_get_N_steps(const char *filename,
				const char *field_name, double *field,
				uint32_t *N_steps)
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
		int readed_step;
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
		      uint32_t N_elems, int *readed_step)
{
	int status = nb_cfreader_check_token(cfr, "Step");
	if (0 != status) {
		if (NB_CFREADER_EOF == status)
			status = STEP_READ_FINISH;
		else
			status = 1;
		goto EXIT;
	}
	
	status = nb_cfreader_read_int(cfr, readed_step);
	if (0 != status)
		goto EXIT;
	
	double time;
	status = nb_cfreader_read_var_double(cfr, "Time", &time);
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
EXIT:
	return status;
}

static int check_field(nb_cfreader_t *cfr, uint32_t N_nodes,
		       uint32_t N_edges, uint32_t N_elems)
{
	int status = nb_cfreader_check_line(cfr, "Field");
	if (0 != status)
		goto EXIT;
	
	char var[30];
	status = nb_cfreader_read_var_token(cfr, "Type", var);
	if (0 != status)
		goto EXIT;

	int N_comp = get_N_components(var);

	if (0 == N_comp) {
		status = 1;
		goto EXIT;
	}

	status = nb_cfreader_read_var_string(cfr, "Name", var);
	if (0 != status)
		goto EXIT;

	status = nb_cfreader_read_var_token(cfr, "Support", var);
	if (0 != status)
		goto EXIT;

	int N = get_field_length(var, N_nodes, N_edges, N_elems);
	if (0 == N) {
		status = 1;
		goto EXIT;
	}

	status = nb_cfreader_check_line(cfr, "Data");
	if (0 != status)
		goto EXIT;
	
	for (uint32_t i = 0; i < N; i++) {
		status = nb_cfreader_jump_line(cfr);
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

int nb_mesh2D_field_read_last(const char *dir_saved_results,
			      const char *field_name, double *field)
{
	/* TEMPORAL */
	return 0;
}
