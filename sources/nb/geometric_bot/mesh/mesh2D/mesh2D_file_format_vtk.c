#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "nb/geometric_bot.h"
#include "nb/cfreader_bot.h"

#include "mesh2D_private.h"
#include "elements2D/msh3trg_vtk.h"
#include "elements2D/mshquad_vtk.h"
#include "elements2D/mshpoly_vtk.h"
#include "elements2D/mshpack_vtk.h"

static void vtk_write_header(FILE *fp, const nb_mesh2D_t *mesh,
			     const char *extra_file);
static void vtk_write_data(FILE *fp, const nb_mesh2D_t *mesh);
static int vtk_get_cell_type(int N_adj);
static int vtk_read_header(nb_cfreader_t *cfr, nb_mesh2D_type *type);
static int get_mesh_type(char var[100], nb_mesh2D_type *type);
static int vtk_read_data(nb_cfreader_t *cfr, nb_mesh2D_t *mesh);
static int vtk_check_data_cell_types(nb_cfreader_t *cfr,
				     const nb_mesh2D_t *mesh);

int nb_mesh2D_save_vtk(const nb_mesh2D_t *mesh, const char *name)
{
	int status = 1;
	FILE *fp = fopen(name, "w");
	if (NULL == fp)
		goto EXIT;

	char extra_file[100];
	sprintf(extra_file, "%s", name);
	char *pch = strstr(extra_file, ".vtk");
	if (NULL == pch)
		goto EXIT;
	sprintf(pch, "_extra.txt");

	vtk_write_header(fp, mesh, extra_file);
	vtk_write_data(fp, mesh);

	fclose(fp);
	status = 0;
EXIT:
	return status;
}

static void vtk_write_header(FILE *fp, const nb_mesh2D_t *mesh,
			     const char *extra_file)
{
	fprintf(fp, "# vtk DataFile Version 2.0\n");
	const char *type = nb_mesh2D_get_type_string(mesh);
	fprintf(fp, "# nbots nb_mesh2D_t 1.0 type=%s ", type);
	fprintf(fp, "extra_file=%s\n", extra_file);
	fprintf(fp, "ASCII\n");
	fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
}

static void vtk_write_data(FILE *fp, const nb_mesh2D_t *mesh)
{
	nb_mesh2D_private_i priv;
	nb_mesh2D_init_private_interface(&priv, mesh);

	uint32_t N_nodes = nb_mesh2D_get_N_nodes(mesh);
	fprintf(fp, "POINTS %i float\n", N_nodes);
	for (uint32_t i = 0; i < N_nodes; i++) {
		float x = nb_mesh2D_node_get_x(mesh, i);
		float y = nb_mesh2D_node_get_y(mesh, i);
		fprintf(fp, " %f %f 0\n", x, y);
	}
	
	uint32_t N_total_adj = priv.get_N_total_adj(mesh->msh);
	uint32_t N_elems = nb_mesh2D_get_N_elems(mesh);
	fprintf(fp, "CELLS %i %i\n", N_elems, N_total_adj + N_elems);
	for (uint32_t i = 0; i < N_elems; i++) {
		int N_adj = nb_mesh2D_elem_get_N_adj(mesh, i);
		fprintf(fp, " %i ", N_adj);
		for (int j = 0; j < N_adj; j++) {
			uint32_t id = nb_mesh2D_elem_get_adj(mesh, i, j);
			fprintf(fp, "%i ", id);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "CELL_TYPES %i\n", N_elems);
	for (uint32_t i = 0; i < N_elems; i++) {
		int N_adj = nb_mesh2D_elem_get_N_adj(mesh, i);
		int type = vtk_get_cell_type(N_adj);
		fprintf(fp, " %i\n", type);
	}
}

static int vtk_get_cell_type(int N_adj)
{
	int type;
	switch (N_adj) {
	case 3:
		type = 5; /* VTK_TRIANGLE */
		break;
	case 4:
		type = 9; /* VTK_QUAD */
		break;
	default:
		type = 7; /* VTK_POLYGON */
	}
	return type;
}

int nb_mesh2D_read_vtk(nb_mesh2D_t *mesh, const char *name)
{
	nb_cfreader_t *cfr = nb_cfreader_create();
	int status = nb_cfreader_open_file(cfr, name);
	if (status != 0)
		goto EXIT;

	nb_mesh2D_type type;
	status = vtk_read_header(cfr, &type);
	if (status != 0)
		goto CLOSE_FILE;

	if (nb_mesh2D_get_type(mesh) != type) {
		status = 1;
		goto CLOSE_FILE;
	}

	status = vtk_read_data(cfr, mesh);
	if (status != 0)
		goto CLOSE_FILE;

	status = 0;
CLOSE_FILE:
	nb_cfreader_close_file(cfr);
EXIT:
	nb_cfreader_destroy(cfr);
	return status;
}

static int vtk_read_header(nb_cfreader_t *cfr, nb_mesh2D_type *type)
{
	int status = 1;
	if (!nb_cfreader_check_line(cfr, "# vtk DataFile Version 2.0"))
		goto EXIT;
	if (!nb_cfreader_check_token(cfr, "# nbots nb_mesh2D_t 1.0"))
		goto EXIT;
	char var[100];
	if (0 != nb_cfreader_read_token(cfr, var))
		goto EXIT;
	if (0 != get_mesh_type(var, type))
		goto EXIT;
	if (!nb_cfreader_check_line(cfr, "ASCII"))
		goto EXIT;
	if (!nb_cfreader_check_line(cfr, "DATASET UNSTRUCTURED_GRID"))
		goto EXIT;
	status = 0;
EXIT:
	return status;
}

static int get_mesh_type(char var[100], nb_mesh2D_type *type)
{
	int status = 0;	
	if (0 == strncmp(var, "NB_TRIAN", 8)) {
		*type = NB_TRIAN;
		goto EXIT;
	}
	if (0 == strncmp(var, "NB_QUAD", 7)) {
		*type = NB_QUAD;
		goto EXIT;
	}
	if (0 == strncmp(var, "NB_POLY", 7)) {
		*type = NB_POLY;
		goto EXIT;
	}
	if (0 == strncmp(var, "NB_DISK", 7)) {
		*type = NB_DISK;
		goto EXIT;
	}
	status = 1;
EXIT:
	return status;
}

static int vtk_read_data(nb_cfreader_t *cfr, nb_mesh2D_t *mesh)
{
	int status;
	switch (nb_mesh2D_get_type(mesh)) {
	case NB_TRIAN:
		status = nb_msh3trg_read_vtk_data(cfr, mesh->msh);
		break;
	case NB_QUAD:
		status = nb_mshquad_read_vtk_data(cfr, mesh->msh);
		break;
	case NB_POLY:
		status = nb_mshpoly_read_vtk_data(cfr, mesh->msh);
		break;
	case NB_DISK:
		status = nb_mshpack_read_vtk_data(cfr, mesh->msh);
		break;
	}
	if (0 != status)
		goto EXIT;
	status = vtk_check_data_cell_types(cfr, mesh);
EXIT:
	return status;
}

static int vtk_check_data_cell_types(nb_cfreader_t *cfr,
				    const nb_mesh2D_t *mesh)
{
	int status = 1;

	char var[100];
	if (0 != nb_cfreader_read_token(cfr, var))
		goto EXIT;

	if (strncmp(var, "CELL_TYPES", 10) != 0)
		goto EXIT;
	
	int N_elems;
	if (0 != nb_cfreader_read_int(cfr, &N_elems))
		goto EXIT;

	if (nb_mesh2D_get_N_elems(mesh) != N_elems)
		goto EXIT;

	for (uint32_t i = 0; i < N_elems; i++) {
		int N_adj = nb_mesh2D_elem_get_N_adj(mesh, i);
		int type = vtk_get_cell_type(N_adj);
		int readed_type;
		if (0 != nb_cfreader_read_int(cfr, &readed_type))
			goto EXIT;
		if (type != readed_type)
			goto EXIT;
	}

	status = 0;
EXIT:
	return status;
}
