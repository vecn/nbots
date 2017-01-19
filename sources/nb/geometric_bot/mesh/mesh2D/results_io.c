#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>

#include "nb/memory_bot.h"
#include "nb/io_bot.h"
#include "nb/geometric_bot.h"

#define STEP_READ_FINISH 100

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
static void draw_nodal_damage(const char *dir_output,
			      const nb_mesh2D_t *mesh,
			      const double *damage,
			      int step);
static void draw_face_damage(const char *dir_output,
			     const nb_mesh2D_t *mesh,
			     const double *disp,
			     const double *damage,
			     int step);

int nb_mesh2D_read_and_draw_results(const char *dir_saved_results,
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

		draw_nodal_damage(dir_output, mesh, damage, step);
		draw_face_damage(dir_output, mesh, disp, damage, step);
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

static void draw_nodal_damage(const char *dir_output,
			      const nb_mesh2D_t *mesh,
			      const double *damage,
			      int step)
{
	uint32_t N_nodes = nb_mesh2D_get_N_nodes(mesh);

	double *node_damage = nb_allocate_mem(N_nodes * sizeof(*node_damage));
	nb_mesh2D_extrapolate_elems_to_nodes(mesh, 1, damage, node_damage);
	
	char name[100];
	sprintf(name, "%s/damage_field_%i.png", dir_output, step);

	nb_mesh2D_export_draw(mesh, name, 1000, 700, NB_NODE, NB_FIELD,
			      node_damage, true);
	nb_free_mem(node_damage);
}

static void draw_face_damage(const char *dir_output,
			     const nb_mesh2D_t *mesh,
			     const double *disp,
			     const double *damage,
			     int step)
{
	/*
	compute_damage(damage, faces, mesh, elem_damage, &glq, &eval_dmg);

	uint32_t N_faces = nb_mesh2D_get_N_edges(mesh);
	double max_dmg = 0;
	for (uint32_t i = 0; i < N_faces; i++) {
		if (max_dmg < face_damage[i])
			max_dmg = face_damage[i];
	}
	for (uint32_t i = 0; i < N_faces; i++)
		face_damage[i] /= max_dmg;
	nb_mesh2D_export_draw(mesh, "./CVFA_dmg.png", 1000, 800,
			      NB_FACE, NB_FIELD,
			      face_damage, true);
	*/
}
