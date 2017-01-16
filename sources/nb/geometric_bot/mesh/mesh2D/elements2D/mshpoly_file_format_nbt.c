#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>

#include "nb/memory_bot.h"
#include "nb/cfreader_bot.h"
#include "nb/geometric_bot.h"

#include "mshpoly_struct.h"
#include "mshpoly_file_format_nbt.h"

static void write_header(FILE *fp, const nb_mshpoly_t *msh);
static void write_nodes(FILE *fp, const nb_mshpoly_t *msh);
static void write_edges(FILE *fp, const nb_mshpoly_t *msh);
static void write_elems(FILE *fp, const nb_mshpoly_t *msh);
static void write_invtx(FILE *fp, const nb_mshpoly_t *msh);
static void write_insgm(FILE *fp, const nb_mshpoly_t *msh);

static int read_header(nb_cfreader_t *cfr, nb_mshpoly_t *msh);
static int read_nodes(nb_cfreader_t *cfr, nb_mshpoly_t *msh);
static int read_edges(nb_cfreader_t *cfr, nb_mshpoly_t *msh);
static int read_elems(nb_cfreader_t *cfr, nb_mshpoly_t *msh);
static int read_invtx(nb_cfreader_t *cfr, nb_mshpoly_t *msh);
static int read_insgm(nb_cfreader_t *cfr, nb_mshpoly_t *msh);

void nb_mshpoly_write_data_nbt(FILE *fp, const void *msh)
{
	write_header(fp, msh);
	write_nodes(fp, msh);
	write_edges(fp, msh);
	write_elems(fp, msh);
	write_invtx(fp, msh);
	write_insgm(fp, msh);
}

static void write_header(FILE *fp, const nb_mshpoly_t *msh)
{
	fprintf(fp, "N_nodes = %i\n", msh->N_nod);
	fprintf(fp, "N_edges = %i\n", msh->N_edg);
	fprintf(fp, "N_elems = %i\n", msh->N_elems);
	fprintf(fp, "N_input_vtx = %i\n", msh->N_vtx);
	fprintf(fp, "N_input_sgm = %i\n", msh->N_sgm);
}

static void write_nodes(FILE *fp, const nb_mshpoly_t *msh)
{
	fprintf(fp, "Nodes\n");
	for (uint32_t i = 0; i < msh->N_nod; i++)
		fprintf(fp, "  %e  %e\n", msh->nod[i * 2], msh->nod[i*2+1]);
	fprintf(fp, "End Nodes\n");
}

static void write_edges(FILE *fp, const nb_mshpoly_t *msh)
{
	fprintf(fp, "Edges\n");
	for (uint32_t i = 0; i < msh->N_edg; i++)
		fprintf(fp, "  %i %i\n",
			msh->edg[i * 2] + 1,
			msh->edg[i*2+1] + 1);
	fprintf(fp, "End Edges\n");
}

static void write_elems(FILE *fp, const nb_mshpoly_t *msh)
{
	fprintf(fp, "Elems\n");
	fprintf(fp, "  Centroids\n");
	for (uint32_t i = 0; i < msh->N_elems; i++)
		fprintf(fp, "    %e %e\n", msh->cen[i * 2], msh->cen[i*2+1]);
	fprintf(fp, "  End Centroids\n");

	fprintf(fp, "  Types\n");
	for (uint32_t i = 0; i < msh->N_elems; i++)
		fprintf(fp, "    %i\n ", msh->N_adj[i]);
	fprintf(fp, "  End Types\n");

	fprintf(fp, "  Adjacencies\n");
	for (uint32_t i = 0; i < msh->N_elems; i++) {
		for (uint32_t j = 0; j < msh->N_adj[i]; j++)
			fprintf(fp, "%i ", msh->adj[i][j] + 1);
		fprintf(fp, "\n");
	}
	fprintf(fp, "  End Adjacencies\n");

	fprintf(fp, "  Neighbors\n");
	fprintf(fp, "  # The null neighbors at the boundary have an ID\n");
	fprintf(fp, "  # greater than the number of elements.\n");
	for (uint32_t i = 0; i < msh->N_elems; i++) {
		for (uint32_t j = 0; j < msh->N_adj[i]; j++)
			fprintf(fp, "%i ", msh->ngb[i][j] + 1);
		fprintf(fp, "\n");
	}
	fprintf(fp, "  End Neighbours\n");

	fprintf(fp, "End Elems\n");
}

static void write_invtx(FILE *fp, const nb_mshpoly_t *msh)
{
	fprintf(fp, "Input vertices\n");
	for (uint32_t i = 0; i < msh->N_vtx; i++)
		fprintf(fp, " %i ->  %i\n", i + 1, msh->vtx[i] + 1);
	fprintf(fp, "End input vertices\n");
}

static void write_insgm(FILE *fp, const nb_mshpoly_t *msh)
{
	fprintf(fp, "Input segments\n");
	for (uint32_t i = 0; i < msh->N_sgm; i++) {
		fprintf(fp, "    %i", msh->N_nod_x_sgm[i]);
		for (uint32_t j = 0; j < msh->N_nod_x_sgm[i]; j++)
			fprintf(fp, "%i ", msh->nod_x_sgm[i][j] + 1);
		fprintf(fp, "\n");
	}
	fprintf(fp, "End input segments\n");
}

int nb_mshpoly_read_data_nbt(nb_cfreader_t *cfr, void *msh_ptr)
{
	nb_mshpoly_t *msh = msh_ptr;

	int status = read_header(cfr, msh);
	if (0 != status)
		goto EXIT;

	nb_mshpoly_set_arrays_memory(msh);

	status = read_nodes(cfr, msh);
	if (0 != status)
		goto CLEAN_MEMBLOCK;

	status = read_edges(cfr, msh);
	if (0 != status)
		goto CLEAN_MEMBLOCK;

	status = read_elems(cfr, msh);
	if (0 != status)
		goto CLEAN_MEMBLOCK;

	status = read_invtx(cfr, msh);
	if (0 != status)
		goto CLEAN_MEMBLOCK;

	status = read_insgm(cfr, msh);
	if (0 != status)
		goto CLEAN_MEMBLOCK;
EXIT:
	return status;
CLEAN_MEMBLOCK:
	nb_free_mem(msh->nod);
	return status;
}

static int read_header(nb_cfreader_t *cfr, nb_mshpoly_t *msh)
{
	uint32_t N;
	int status = nb_cfreader_read_var_uint(cfr, "N_nodes", &N);
	if (0 != status)
		goto EXIT;
	msh->N_nod = N;

	status = nb_cfreader_read_var_uint(cfr, "N_edges", &N);
	if (0 != status)
		goto EXIT;
	N = msh->N_edg;

	status = nb_cfreader_read_var_uint(cfr, "N_elems", &N);
	if (0 != status)
		goto EXIT;
	N = msh->N_elems;

	status = nb_cfreader_read_var_uint(cfr, "N_input_vtx", &N);
	if (0 != status)
		goto EXIT;
	N = msh->N_vtx;

	status = nb_cfreader_read_var_uint(cfr, "N_input_sgm", &N);
	if (0 != status)
		goto EXIT;
	N = msh->N_sgm;

EXIT:
	return status;
}

static int read_nodes(nb_cfreader_t *cfr, nb_mshpoly_t *msh)
{
	int status = nb_cfreader_check_line(cfr, "Nodes");
	if (0 != status)
		goto EXIT;

	for (uint32_t i = 0; i < msh->N_nod; i++) {
		double val;
		status = nb_cfreader_read_double(cfr, &val);
		if (0 != status)
			goto EXIT;
		msh->nod[i * 2] = val;

		status = nb_cfreader_read_double(cfr, &val);
		if (0 != status)
			goto EXIT;
		msh->nod[i*2+1] = val;
	}

	status = nb_cfreader_check_line(cfr, "End Nodes");
	if (0 != status)
		goto EXIT;
EXIT:
	return status;
}

static int read_edges(nb_cfreader_t *cfr, nb_mshpoly_t *msh)
{
	int status = nb_cfreader_check_line(cfr, "Edges");
	if (0 != status)
		goto EXIT;

	for (uint32_t i = 0; i < msh->N_edg; i++) {
		uint32_t val;
		status = nb_cfreader_read_uint(cfr, &val);
		if (0 != status)
			goto EXIT;
		msh->edg[i * 2] = val - 1;

		status = nb_cfreader_read_uint(cfr, &val);
		if (0 != status)
			goto EXIT;
		msh->edg[i*2+1] = val - 1;
	}

	status = nb_cfreader_check_line(cfr, "End Edges");
	if (0 != status)
		goto EXIT;
EXIT:
	return status;
}

static int read_elems(nb_cfreader_t *cfr, nb_mshpoly_t *msh)
{
	return 1;
}

static int read_invtx(nb_cfreader_t *cfr, nb_mshpoly_t *msh)
{
	return 1;
}

static int read_insgm(nb_cfreader_t *cfr, nb_mshpoly_t *msh)
{
	return 1;
}
