#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>

#include "nb/memory_bot.h"
#include "nb/io_bot.h"
#include "nb/geometric_bot.h"

#include "mshpoly_struct.h"
#include "mshpoly_file_format_nbt.h"

#define NNODES_VARNAME "N_nodes"
#define NEDGES_VARNAME "N_edges"
#define NELEMS_VARNAME "N_elems"
#define NINVTX_VARNAME "N_input_vtx"
#define NINSGM_VARNAME "N_input_sgm"

#define NODES_OPEN_SECTION "Nodes"
#define NODES_CLOSE_SECTION "End Nodes"

#define EDGES_OPEN_SECTION "Edges"
#define EDGES_CLOSE_SECTION "End Edges"

#define ELEMS_OPEN_SECTION "Elems"
#define ELEMS_CLOSE_SECTION "End Elems"

#define ELEMS_CEN_OPEN_SUBSECTION "Centroids"
#define ELEMS_CEN_CLOSE_SUBSECTION "End Centroids"

#define ELEMS_NSIDES_OPEN_SUBSECTION "N-sides"
#define ELEMS_NSIDES_CLOSE_SUBSECTION "End N-sides"

#define ELEMS_ADJ_OPEN_SUBSECTION "Adjacencies"
#define ELEMS_ADJ_CLOSE_SUBSECTION "End Adjacencies"

#define ELEMS_NGB_OPEN_SUBSECTION "Neighbors"
#define ELEMS_NGB_CLOSE_SUBSECTION "End Neighbors"

#define INVTX_OPEN_SECTION "Input Vertices"
#define INVTX_CLOSE_SECTION "End Input Vertices"

#define INSGM_OPEN_SECTION "Input Segments"
#define INSGM_CLOSE_SECTION "End Input Segments"

#define INSGM_SIZE_OPEN_SUBSECTION "N-nodes"
#define INSGM_SIZE_CLOSE_SUBSECTION "End N-nodes"

#define INSGM_IDNODE_OPEN_SUBSECTION "ID-nodes"
#define INSGM_IDNODE_CLOSE_SUBSECTION "End ID-nodes"

static void write_header(FILE *fp, const nb_mshpoly_t *msh);
static void write_nodes(FILE *fp, const nb_mshpoly_t *msh);
static void write_edges(FILE *fp, const nb_mshpoly_t *msh);
static void write_elems(FILE *fp, const nb_mshpoly_t *msh);
static void write_elems_centroids(FILE *fp, const nb_mshpoly_t *msh);
static void write_elems_N_sides(FILE *fp, const nb_mshpoly_t *msh);
static void write_elems_adj(FILE *fp, const nb_mshpoly_t *msh);
static void write_elems_ngb(FILE *fp, const nb_mshpoly_t *msh);
static void write_invtx(FILE *fp, const nb_mshpoly_t *msh);
static void write_insgm(FILE *fp, const nb_mshpoly_t *msh);
static void write_insgm_size(FILE *fp, const nb_mshpoly_t *msh);
static void write_insgm_idvtx(FILE *fp, const nb_mshpoly_t *msh);

static int read_header(nb_cfreader_t *cfr, nb_mshpoly_t *msh);
static int read_nodes(nb_cfreader_t *cfr, nb_mshpoly_t *msh);
static int read_edges(nb_cfreader_t *cfr, nb_mshpoly_t *msh);
static int read_elems(nb_cfreader_t *cfr, nb_mshpoly_t *msh);
static int read_elems_centroids(nb_cfreader_t *cfr, nb_mshpoly_t *msh);
static int read_elems_N_sides(nb_cfreader_t *cfr, nb_mshpoly_t *msh);
static int read_elems_adj(nb_cfreader_t *cfr, nb_mshpoly_t *msh);
static int read_elems_adj_id(nb_cfreader_t *cfr, nb_mshpoly_t *msh,
			     uint32_t i, int j);
static int read_elems_ngb(nb_cfreader_t *cfr, nb_mshpoly_t *msh);
static int read_elems_ngb_id(nb_cfreader_t *cfr, nb_mshpoly_t *msh,
			     uint32_t i, int j);
static int read_invtx(nb_cfreader_t *cfr, nb_mshpoly_t *msh);
static int read_invtx_node(nb_cfreader_t *cfr, nb_mshpoly_t *msh, uint32_t i);
static int read_insgm(nb_cfreader_t *cfr, nb_mshpoly_t *msh);
static int read_insgm_size(nb_cfreader_t *cfr, nb_mshpoly_t *msh);
static int read_insgm_idvtx(nb_cfreader_t *cfr, nb_mshpoly_t *msh);
static int read_insgm_idvtx_node(nb_cfreader_t *cfr, nb_mshpoly_t *msh,
				 uint32_t i, uint32_t j);

void nb_mshpoly_write_data_nbt(FILE *fp, const void *msh)
{
	write_header(fp, msh);
	fprintf(fp, "\n");

	write_nodes(fp, msh);
	fprintf(fp, "\n");

	write_edges(fp, msh);
	fprintf(fp, "\n");

	write_elems(fp, msh);
	fprintf(fp, "\n");

	write_invtx(fp, msh);
	fprintf(fp, "\n");

	write_insgm(fp, msh);
}

static void write_header(FILE *fp, const nb_mshpoly_t *msh)
{
	fprintf(fp, "%s = %i\n", NNODES_VARNAME, msh->N_nod);
	fprintf(fp, "%s = %i\n", NEDGES_VARNAME, msh->N_edg);
	fprintf(fp, "%s = %i\n", NELEMS_VARNAME, msh->N_elems);
	fprintf(fp, "%s = %i\n", NINVTX_VARNAME, msh->N_vtx);
	fprintf(fp, "%s = %i\n", NINSGM_VARNAME, msh->N_sgm);
}

static void write_nodes(FILE *fp, const nb_mshpoly_t *msh)
{
  fprintf(fp, "%s\n", NODES_OPEN_SECTION);
	for (uint32_t i = 0; i < msh->N_nod; i++)
		fprintf(fp, "  %e  %e\n", msh->nod[i * 2], msh->nod[i*2+1]);
	fprintf(fp, "%s\n", NODES_CLOSE_SECTION);
}

static void write_edges(FILE *fp, const nb_mshpoly_t *msh)
{
  fprintf(fp, "%s\n", EDGES_OPEN_SECTION);
	for (uint32_t i = 0; i < msh->N_edg; i++)
		fprintf(fp, "  %i %i\n",
			msh->edg[i * 2] + 1,
			msh->edg[i*2+1] + 1);
	fprintf(fp, "%s\n", EDGES_CLOSE_SECTION);
}

static void write_elems(FILE *fp, const nb_mshpoly_t *msh)
{
  fprintf(fp, "%s\n", ELEMS_OPEN_SECTION);

	write_elems_centroids(fp, msh);
	fprintf(fp, "\n");

	write_elems_N_sides(fp, msh);
	fprintf(fp, "\n");

	write_elems_adj(fp, msh);
	fprintf(fp, "\n");

	write_elems_ngb(fp, msh);

	fprintf(fp, "%s\n", ELEMS_CLOSE_SECTION);
}

static void write_elems_centroids(FILE *fp, const nb_mshpoly_t *msh)
{
	fprintf(fp, "  %s\n", ELEMS_CEN_OPEN_SUBSECTION);
	for (uint32_t i = 0; i < msh->N_elems; i++)
		fprintf(fp, "    %e %e\n", msh->cen[i * 2], msh->cen[i*2+1]);
	fprintf(fp, "  %s\n", ELEMS_CEN_CLOSE_SUBSECTION);
}

static void write_elems_N_sides(FILE *fp, const nb_mshpoly_t *msh)
{
	fprintf(fp, "  %s\n", ELEMS_NSIDES_OPEN_SUBSECTION);
	for (uint32_t i = 0; i < msh->N_elems; i++) {
		if (0 == i % 10)
			fprintf(fp, "    ");
		fprintf(fp, "%i  ", msh->N_adj[i]);
		
		if (0 == (i+1) % 10 ||  (i+1) == msh->N_elems)
			fprintf(fp, "\n");
	}
	fprintf(fp, "  %s\n", ELEMS_NSIDES_CLOSE_SUBSECTION);
}

static void write_elems_adj(FILE *fp, const nb_mshpoly_t *msh)
{
	fprintf(fp, "  %s\n", ELEMS_ADJ_OPEN_SUBSECTION);
	for (uint32_t i = 0; i < msh->N_elems; i++) {
		fprintf(fp, "    ");
		for (uint32_t j = 0; j < msh->N_adj[i]; j++)
			fprintf(fp, "%i ", msh->adj[i][j] + 1);
		fprintf(fp, "\n");
	}
	fprintf(fp, "  %s\n", ELEMS_ADJ_CLOSE_SUBSECTION);
}

static void write_elems_ngb(FILE *fp, const nb_mshpoly_t *msh)
{
	fprintf(fp, "  %s\n", ELEMS_NGB_OPEN_SUBSECTION);
	fprintf(fp, "  # ID = 0 indicates null neighbors (boundaries).\n");
	for (uint32_t i = 0; i < msh->N_elems; i++) {
		fprintf(fp, "    ");
		for (uint32_t j = 0; j < msh->N_adj[i]; j++) {
			if (msh->ngb[i][j] < msh->N_elems)
				fprintf(fp, "%i ", msh->ngb[i][j] + 1);
			else
				fprintf(fp, "0 ");
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "  %s\n", ELEMS_NGB_CLOSE_SUBSECTION);
}

static void write_invtx(FILE *fp, const nb_mshpoly_t *msh)
{
	int N_null = 0;
	for (uint32_t i = 0; i < msh->N_vtx; i++) {
		if (msh->vtx[i] < msh->N_vtx)
			N_null += 1;
	}
	if (N_null < msh->N_vtx) {
		fprintf(fp, "%s\n", INVTX_OPEN_SECTION);
		if (N_null > 0)
			fprintf(fp, "  # ID = 0 indicates that mesh" \
				"vertex does not exist.\n");
		for (uint32_t i = 0; i < msh->N_vtx; i++) {
			if (msh->vtx[i] < msh->N_nod)
				fprintf(fp, "  %i <- %i\n", i + 1,
					msh->vtx[i] + 1);
			else
				fprintf(fp, "  %i <- 0\n", i + 1);
		}
		fprintf(fp, "%s\n", INVTX_CLOSE_SECTION);
	}
}

static void write_insgm(FILE *fp, const nb_mshpoly_t *msh)
{
	int N_null = 0;
	for (uint32_t i = 0; i < msh->N_sgm; i++) {
		if (0 == msh->N_nod_x_sgm[i])
			N_null += 1;
	}
	if (N_null < msh->N_sgm) {
		fprintf(fp, "%s\n", INSGM_OPEN_SECTION);
		if (N_null > 0)
			fprintf(fp, "  # Size = 0 indicates that mesh" \
				"segment does not exist.\n");
		write_insgm_size(fp, msh);
		fprintf(fp, "\n");
		write_insgm_idvtx(fp, msh);
		fprintf(fp, "%s\n", INSGM_CLOSE_SECTION);
	}
}

static void write_insgm_size(FILE *fp, const nb_mshpoly_t *msh)
{
	fprintf(fp, "  %s\n", INSGM_SIZE_OPEN_SUBSECTION);
	for (uint32_t i = 0; i < msh->N_sgm; i++) {
		if (0 == i % 10)
			fprintf(fp, "    ");
		fprintf(fp, "%i ", msh->N_nod_x_sgm[i]);
		
		if (0 == (i+1) % 10 ||  (i+1) == msh->N_elems)
			fprintf(fp, "\n");
	}
	fprintf(fp, "  %s\n", INSGM_SIZE_CLOSE_SUBSECTION);
}

static void write_insgm_idvtx(FILE *fp, const nb_mshpoly_t *msh)
{
	fprintf(fp, "  %s\n", INSGM_IDNODE_OPEN_SUBSECTION);
	for (uint32_t i = 0; i < msh->N_sgm; i++) {
		if (msh->N_nod_x_sgm[i] > 0) {
			fprintf(fp, "   [ ");
			for (uint32_t j = 0; j < msh->N_nod_x_sgm[i]; j++) {
				if (0 == j % 10 && j > 0)
					fprintf(fp, "     ");
				fprintf(fp, "%i ", msh->nod_x_sgm[i][j] + 1);
		
				if ((j+1) == msh->N_nod_x_sgm[i]) {
					fprintf(fp, "]\n");
				} else {
					fprintf(fp, "-> ");					
					if (0 == (j+1) % 10)
						fprintf(fp, "\n");
				}
			}
		} else {
			fprintf(fp, "    []\n");
		}
		if (i < msh->N_sgm - 1)
			fprintf(fp, "\n");
	}
	fprintf(fp, "  %s\n", INSGM_IDNODE_CLOSE_SUBSECTION);
}

int nb_mshpoly_read_data_nbt(nb_cfreader_t *cfr, void *msh_ptr)
{
	nb_mshpoly_t *msh = msh_ptr;

	msh->nod = NULL;

	int status = read_header(cfr, msh);
	if (0 != status)
		goto EXIT;

	nb_mshpoly_set_arrays_memory(msh);

	msh->adj[0] = NULL;
	msh->nod_x_sgm[0] = NULL;

	status = read_nodes(cfr, msh);
	if (0 != status)
		goto EXIT;

	status = read_edges(cfr, msh);
	if (0 != status)
		goto EXIT;

	status = read_elems(cfr, msh);
	if (0 != status)
		goto EXIT;

	status = read_invtx(cfr, msh);
	if (0 != status)
		goto EXIT;

	status = read_insgm(cfr, msh);
	if (0 != status)
		goto EXIT;

	status = NB_MSHPOLY_READ_SUCCESS;
EXIT:
	if (0 != status) {
		if (NULL != msh->nod) {
			if (NULL != msh->adj[0])
				nb_free_mem(msh->adj[0]);
			if (NULL != msh->nod_x_sgm[0])
			nb_free_mem(msh->nod_x_sgm[0]);
			nb_free_mem(msh->nod);
		}
	}
	return status;
}

static int read_header(nb_cfreader_t *cfr, nb_mshpoly_t *msh)
{
	uint32_t N;
	int status = nb_cfreader_read_var_uint(cfr, NNODES_VARNAME, &N);
	if (0 != status) {
		status = NB_MSHPOLY_READ_ERROR_NNODES_VARNAME;
		goto EXIT;
	}
	msh->N_nod = N;

	status = nb_cfreader_read_var_uint(cfr, NEDGES_VARNAME, &N);
	if (0 != status) {
		status = NB_MSHPOLY_READ_ERROR_NEDGES_VARNAME;
		goto EXIT;
	}
	msh->N_edg = N;

	status = nb_cfreader_read_var_uint(cfr, NELEMS_VARNAME, &N);
	if (0 != status) {
		status = NB_MSHPOLY_READ_ERROR_NELEMS_VARNAME;
		goto EXIT;
	}
	msh->N_elems = N;

	status = nb_cfreader_read_var_uint(cfr, NINVTX_VARNAME, &N);
	if (0 != status) {
		status = NB_MSHPOLY_READ_ERROR_NINVTX_VARNAME;
		goto EXIT;
	}
	msh->N_vtx = N;

	status = nb_cfreader_read_var_uint(cfr, NINSGM_VARNAME, &N);
	if (0 != status) {
		status = NB_MSHPOLY_READ_ERROR_NINSGM_VARNAME;
		goto EXIT;
	}
	msh->N_sgm = N;

EXIT:
	return status;
}

static int read_nodes(nb_cfreader_t *cfr, nb_mshpoly_t *msh)
{
	int status = nb_cfreader_check_line(cfr, NODES_OPEN_SECTION);
	if (0 != status) {
		status = NB_MSHPOLY_READ_ERROR_NODES_OPEN_SECTION;
		goto EXIT;
	}

	for (uint32_t i = 0; i < msh->N_nod; i++) {
		double val;
		status = nb_cfreader_read_double(cfr, &val);
		if (0 != status) {
			status = NB_MSHPOLY_READ_ERROR_NODES_BAD_INPUT;
			goto EXIT;
		}
		msh->nod[i * 2] = val;

		status = nb_cfreader_read_double(cfr, &val);
		if (0 != status) {
			status = NB_MSHPOLY_READ_ERROR_NODES_BAD_INPUT;
			goto EXIT;
		}
		msh->nod[i*2+1] = val;
	}

	status = nb_cfreader_check_line(cfr, NODES_CLOSE_SECTION);
	if (0 != status) {
		status = NB_MSHPOLY_READ_ERROR_NODES_CLOSE_SECTION;
		goto EXIT;
	}
EXIT:
	return status;
}

static int read_edges(nb_cfreader_t *cfr, nb_mshpoly_t *msh)
{
	int status = nb_cfreader_check_line(cfr, EDGES_OPEN_SECTION);
	if (0 != status) {
		status = NB_MSHPOLY_READ_ERROR_EDGES_OPEN_SECTION;
		goto EXIT;
	}

	for (uint32_t i = 0; i < msh->N_edg; i++) {
		uint32_t val;
		status = nb_cfreader_read_uint(cfr, &val);
		if (0 != status) {
			status = NB_MSHPOLY_READ_ERROR_EDGES_BAD_INPUT;
			goto EXIT;
		}
		if (0 == val || val > msh->N_nod) {
			status = NB_MSHPOLY_READ_ERROR_EDGES_OUT_OF_RANGE;
			goto EXIT;
		}
		msh->edg[i * 2] = val - 1;

		status = nb_cfreader_read_uint(cfr, &val);
		if (0 != status) {
			status = NB_MSHPOLY_READ_ERROR_EDGES_BAD_INPUT;
			goto EXIT;
		}
		if (0 == val || val > msh->N_nod) {
			status = NB_MSHPOLY_READ_ERROR_EDGES_OUT_OF_RANGE;
			goto EXIT;
		}
		msh->edg[i*2+1] = val - 1;
	}

	status = nb_cfreader_check_line(cfr, EDGES_CLOSE_SECTION);
	if (0 != status) {
		status = NB_MSHPOLY_READ_ERROR_EDGES_CLOSE_SECTION;
		goto EXIT;
	}
EXIT:
	return status;
}

static int read_elems(nb_cfreader_t *cfr, nb_mshpoly_t *msh)
{
	int status = nb_cfreader_check_line(cfr, ELEMS_OPEN_SECTION);
	if (0 != status) {
		status = NB_MSHPOLY_READ_ERROR_ELEMS_OPEN_SECTION;
		goto EXIT;
	}
	status = read_elems_centroids(cfr, msh);
	if (0 != status)
		goto EXIT;
	status = read_elems_N_sides(cfr, msh);
	if (0 != status)
		goto EXIT;

	nb_mshpoly_set_mem_of_adj_and_ngb(msh);

	status = read_elems_adj(cfr, msh);
	if (0 != status)
		goto EXIT;
	status = read_elems_ngb(cfr, msh);
	if (0 != status)
		goto EXIT;

	status = nb_cfreader_check_line(cfr, ELEMS_CLOSE_SECTION);
	if (0 != status) {
		status = NB_MSHPOLY_READ_ERROR_ELEMS_CLOSE_SECTION;
		goto EXIT;
	}
EXIT:
	return status;
}

static int read_elems_centroids(nb_cfreader_t *cfr, nb_mshpoly_t *msh)
{
	int status = nb_cfreader_check_line(cfr, ELEMS_CEN_OPEN_SUBSECTION);
	if (0 != status) {
		status = NB_MSHPOLY_READ_ERROR_ELEMS_CEN_OPEN_SUBSECTION;
		goto EXIT;
	}

	for (uint32_t i = 0; i < msh->N_elems; i++) {
		double val;
		status = nb_cfreader_read_double(cfr, &val);
		if (0 != status) {
			status = NB_MSHPOLY_READ_ERROR_ELEMS_CEN_BAD_INPUT;
			goto EXIT;
		}
		msh->cen[i * 2] = val;

		status = nb_cfreader_read_double(cfr, &val);
		if (0 != status) {
			status = NB_MSHPOLY_READ_ERROR_ELEMS_CEN_BAD_INPUT;
			goto EXIT;
		}
		msh->cen[i*2+1] = val;
	}
	status = nb_cfreader_check_line(cfr, ELEMS_CEN_CLOSE_SUBSECTION);
	if (0 != status) {
		status = NB_MSHPOLY_READ_ERROR_ELEMS_CEN_CLOSE_SUBSECTION;
		goto EXIT;
	}
EXIT:
	return status;
}

static int read_elems_N_sides(nb_cfreader_t *cfr, nb_mshpoly_t *msh)
{
	int status = nb_cfreader_check_line(cfr, ELEMS_NSIDES_OPEN_SUBSECTION);
	if (0 != status) {
		status = NB_MSHPOLY_READ_ERROR_ELEMS_NSIDES_OPEN_SUBSECTION;
		goto EXIT;
	}

	for (uint32_t i = 0; i < msh->N_elems; i++) {
		uint32_t val;
		status = nb_cfreader_read_uint(cfr, &val);
		if (0 != status) {
			status = NB_MSHPOLY_READ_ERROR_ELEMS_NSIDES_BAD_INPUT;
			goto EXIT;
		}
		if (0 == val) {
			status = NB_MSHPOLY_READ_ERROR_ELEMS_NSIDES_OUT_RANGE;
			goto EXIT;
		}
		msh->N_adj[i] = val;
	}

	status = nb_cfreader_check_line(cfr, ELEMS_NSIDES_CLOSE_SUBSECTION);
	if (0 != status) {
		status = NB_MSHPOLY_READ_ERROR_ELEMS_NSIDES_CLOSE_SUBSECTION;
		goto EXIT;
	}
EXIT:
	return status;
}

static int read_elems_adj(nb_cfreader_t *cfr, nb_mshpoly_t *msh)
{
	int status = nb_cfreader_check_line(cfr, ELEMS_ADJ_OPEN_SUBSECTION);
	if (0 != status) {
		status = NB_MSHPOLY_READ_ERROR_ELEMS_ADJ_OPEN_SUBSECTION;
		goto EXIT;
	}

	for (uint32_t i = 0; i < msh->N_elems; i++) {
		for (int j = 0; j < msh->N_adj[i]; j++) {
			status = read_elems_adj_id(cfr, msh, i, j);
			if (0 != status)
				goto EXIT;
		}
	}

	status = nb_cfreader_check_line(cfr, ELEMS_ADJ_CLOSE_SUBSECTION);
	if (0 != status) {
		status = NB_MSHPOLY_READ_ERROR_ELEMS_ADJ_CLOSE_SUBSECTION;
		goto EXIT;
	}
EXIT:
	return status;
}

static int read_elems_adj_id(nb_cfreader_t *cfr, nb_mshpoly_t *msh,
			     uint32_t i, int j)
{
	uint32_t val;
	int status = nb_cfreader_read_uint(cfr, &val);
	if (0 != status) {
		status = NB_MSHPOLY_READ_ERROR_ELEMS_ADJ_BAD_INPUT;
		goto EXIT;
	}
	if (0 == val || val > msh->N_nod) {
		status = NB_MSHPOLY_READ_ERROR_ELEMS_NSIDES_OUT_RANGE;
		goto EXIT;
	}
	msh->adj[i][j] = val - 1;
	status = 0;
EXIT:
	return status;
}

static int read_elems_ngb(nb_cfreader_t *cfr, nb_mshpoly_t *msh)
{
	int status = nb_cfreader_check_line(cfr, ELEMS_NGB_OPEN_SUBSECTION);
	if (0 != status) {
		status = NB_MSHPOLY_READ_ERROR_ELEMS_NGB_OPEN_SUBSECTION;
		goto EXIT;
	}

	for (uint32_t i = 0; i < msh->N_elems; i++) {
		for (int j = 0; j < msh->N_adj[i]; j++) {
			status = read_elems_ngb_id(cfr, msh, i, j);
			if (0 != status)
				goto EXIT;
		}
	}

	status = nb_cfreader_check_line(cfr, ELEMS_NGB_CLOSE_SUBSECTION);
	if (0 != status) {
		status = NB_MSHPOLY_READ_ERROR_ELEMS_NGB_CLOSE_SUBSECTION;
		goto EXIT;
	}
EXIT:
	return status;
}

static int read_elems_ngb_id(nb_cfreader_t *cfr, nb_mshpoly_t *msh,
			     uint32_t i, int j)
{
	uint32_t val;
	int status = nb_cfreader_read_uint(cfr, &val);
	if (0 != status) {
		status = NB_MSHPOLY_READ_ERROR_ELEMS_NGB_BAD_INPUT;
		goto EXIT;
	}

	if (0 == val) {
		msh->ngb[i][j] = msh->N_elems;
	} else if (val <= msh->N_elems) {
		msh->ngb[i][j] = val - 1;
	} else {
		status = NB_MSHPOLY_READ_ERROR_ELEMS_NGB_OUT_RANGE;
		goto EXIT;
	}
	status = 0;
EXIT:
	return status;
}

static int read_invtx(nb_cfreader_t *cfr, nb_mshpoly_t *msh)
{
	int status = nb_cfreader_check_line(cfr, INVTX_OPEN_SECTION);
	if (0 != status) {
		if (msh->N_vtx == 0)
			status = 0;
		else
			status = NB_MSHPOLY_READ_ERROR_INVTX_OPEN_SECTION;
		goto EXIT;
	} else {
		if (msh->N_vtx == 0) {
			status = NB_MSHPOLY_READ_ERROR_INVTX_NOT_DEFINED;
			goto EXIT;
		}

	}

	for (uint32_t i = 0; i < msh->N_vtx; i++) {
		status = read_invtx_node(cfr, msh, i);
		if (0 != status)
			goto EXIT;
	}

	status = nb_cfreader_check_line(cfr, INVTX_CLOSE_SECTION);
	if (0 != status) {
		status = NB_MSHPOLY_READ_ERROR_INVTX_CLOSE_SECTION;
		goto EXIT;
	}
EXIT:
	return status;
}

static int read_invtx_node(nb_cfreader_t *cfr, nb_mshpoly_t *msh, uint32_t i)
{
	char val1[15];
	char val2[15];
	int status = nb_cfreader_read_tuple(cfr, val1, val2);
	if (0 != status) {
		status = NB_MSHPOLY_READ_ERROR_INVTX_BAD_INPUT;
		goto EXIT;
	}
		
	char *pch;
	unsigned int N = strtoul(val1, &pch, 10);
	if (pch == val1) {
		status = NB_MSHPOLY_READ_ERROR_INVTX_BAD_INPUT;
		goto EXIT;
	}

	if (0 == N || N > msh->N_vtx) {
		status = NB_MSHPOLY_READ_ERROR_INVTX_OUT_OF_RANGE;
		goto EXIT;
	}
	uint32_t id = N-1;

	N = strtoul(val2, &pch, 10);
	if (pch == val1) {
		status = NB_MSHPOLY_READ_ERROR_INVTX_BAD_INPUT;
		goto EXIT;
	}

	if (N > 0 && N <= msh->N_nod) {
		msh->vtx[id] = N - 1;
	} else if (0 == N) {
		msh->vtx[id] = msh->N_nod;
	} else {
		status = NB_MSHPOLY_READ_ERROR_INVTX_OUT_OF_RANGE;
		goto EXIT;
	}
	status = 0;
EXIT:
	return status;
}

static int read_insgm(nb_cfreader_t *cfr, nb_mshpoly_t *msh)
{
	int status = nb_cfreader_check_line(cfr, INSGM_OPEN_SECTION);
	if (0 != status) {
		if (msh->N_sgm == 0)
			status = 0;
		else
			status = NB_MSHPOLY_READ_ERROR_INSGM_OPEN_SECTION;
		goto EXIT;
	} else {
		if (msh->N_sgm == 0) {
			status = NB_MSHPOLY_READ_ERROR_INSGM_NOT_DEFINED;
			goto EXIT;
		}

	}

	status = read_insgm_size(cfr, msh);
	if (0 != status)
		goto EXIT;

	nb_mshpoly_set_mem_of_nod_x_sgm(msh);

	status = read_insgm_idvtx(cfr, msh);
	if (0 != status)
		goto EXIT;

	status = nb_cfreader_check_line(cfr, INSGM_CLOSE_SECTION);
	if (0 != status) {
		status = NB_MSHPOLY_READ_ERROR_INSGM_CLOSE_SECTION;
		goto EXIT;
	}
EXIT:
	return status;
}

static int read_insgm_size(nb_cfreader_t *cfr, nb_mshpoly_t *msh)
{
	int status = nb_cfreader_check_line(cfr, INSGM_SIZE_OPEN_SUBSECTION);
	if (0 != status) {
		status = NB_MSHPOLY_READ_ERROR_INSGM_SIZE_OPEN_SUBSECTION;
		goto EXIT;
	}

	for (uint32_t i = 0; i < msh->N_sgm; i++) {
		uint32_t val;
		status = nb_cfreader_read_uint(cfr, &val);
		if (0 != status) {
			status = NB_MSHPOLY_READ_ERROR_INSGM_SIZE_BAD_INPUT;
			goto EXIT;
		}

		msh->N_nod_x_sgm[i] = val;
	}

	status = nb_cfreader_check_line(cfr, INSGM_SIZE_CLOSE_SUBSECTION);
	if (0 != status) {
		status = NB_MSHPOLY_READ_ERROR_INSGM_SIZE_CLOSE_SUBSECTION;
		goto EXIT;
	}
EXIT:
	return status;

}

static int read_insgm_idvtx(nb_cfreader_t *cfr, nb_mshpoly_t *msh)
{
	int status = nb_cfreader_check_line(cfr, INSGM_IDNODE_OPEN_SUBSECTION);
	if (0 != status) {
		status = NB_MSHPOLY_READ_ERROR_INSGM_IDS_OPEN_SUBSECTION;
		goto EXIT;
	}

	for (uint32_t i = 0; i < msh->N_sgm; i++) {
		int status = nb_cfreader_check_token(cfr, "[");
		if (0 != status) {
			status = NB_MSHPOLY_READ_ERROR_INSGM_IDS_BAD_OPEN_BRK;
			goto EXIT;
		}
		for (uint32_t j = 0; j < msh->N_nod_x_sgm[i]; j++) {
			status = read_insgm_idvtx_node(cfr, msh, i, j);
			if (0 != status)
				goto EXIT;
		}
		status = nb_cfreader_check_token(cfr, "]");
		if (0 != status) {
			status = NB_MSHPOLY_READ_ERROR_INSGM_IDS_BAD_CLOSE_BRK;
			goto EXIT;
		}
	}

	status = nb_cfreader_check_line(cfr, INSGM_IDNODE_CLOSE_SUBSECTION);
	if (0 != status) {
		status = NB_MSHPOLY_READ_ERROR_INSGM_IDS_CLOSE_SUBSECTION;
		goto EXIT;
	}
EXIT:
	return status;
}

static int read_insgm_idvtx_node(nb_cfreader_t *cfr, nb_mshpoly_t *msh,
				 uint32_t i, uint32_t j)
{
	int status;
	if (j > 0) {
		status = nb_cfreader_check_token(cfr, "->");
		if (0 != status) {
			status = NB_MSHPOLY_READ_ERROR_INSGM_IDS_BAD_LINK;
			goto EXIT;
		}
	}
	uint32_t val;
	status = nb_cfreader_read_uint(cfr, &val);
	if (0 != status) {
		status = NB_MSHPOLY_READ_ERROR_INSGM_IDS_BAD_INPUT;
		goto EXIT;
	}

	if (0 == val || val > msh->N_nod) {
		status = NB_MSHPOLY_READ_ERROR_INSGM_IDS_OUT_OF_RANGE;
		goto EXIT;
	}

	msh->nod_x_sgm[i][j] = val - 1;
	status = 0;
EXIT:
	return status;
}

void nb_mshpoly_get_error_message(int error, char **msg)
{
	switch (error) {
	case NB_MSHPOLY_READ_SUCCESS:
		*msg = "Success";
		break;
	case NB_MSHPOLY_READ_ERROR_NNODES_VARNAME:
		*msg = NNODES_VARNAME		\
			" is not defined";
		break;
	case NB_MSHPOLY_READ_ERROR_NEDGES_VARNAME:
		*msg = NEDGES_VARNAME		\
			" is not defined";
		break;
	case NB_MSHPOLY_READ_ERROR_NELEMS_VARNAME:
		*msg = NELEMS_VARNAME		\
			" is not defined";
		break;
	case NB_MSHPOLY_READ_ERROR_NINVTX_VARNAME:
		*msg = NINVTX_VARNAME		\
			" is not defined";
		break;
	case NB_MSHPOLY_READ_ERROR_NINSGM_VARNAME:
		*msg = NINSGM_VARNAME		\
			" is not defined";
		break;
	case NB_MSHPOLY_READ_ERROR_NODES_OPEN_SECTION:
		*msg = "Nodes section opening is not defined";
		break;
	case NB_MSHPOLY_READ_ERROR_NODES_BAD_INPUT:
		*msg = "Bad coordinates definition in nodes section";
		break;
	case NB_MSHPOLY_READ_ERROR_NODES_CLOSE_SECTION:
		*msg = "Nodes section closing is not defined";
		break;
	case NB_MSHPOLY_READ_ERROR_EDGES_OPEN_SECTION:
		*msg = "Edges section opening is not defined";
		break;
	case NB_MSHPOLY_READ_ERROR_EDGES_BAD_INPUT:
		*msg = "Bad node ID definition in edges section";
		break;
	case NB_MSHPOLY_READ_ERROR_EDGES_OUT_OF_RANGE:
		*msg = "Node ID is out of range in edges section";
		break;
	case NB_MSHPOLY_READ_ERROR_EDGES_CLOSE_SECTION:
		*msg = "Edges section closing is not defined";
		break;
	case NB_MSHPOLY_READ_ERROR_ELEMS_OPEN_SECTION:
		*msg = "Elements section opening is not defined";
		break;
	case NB_MSHPOLY_READ_ERROR_ELEMS_CEN_OPEN_SUBSECTION:
		*msg = "Elements centroids subsection opening is not defined";
		break;
	case NB_MSHPOLY_READ_ERROR_ELEMS_CEN_BAD_INPUT:
		*msg = "Bad coordinates in elements centroids subsection";
		break;
	case NB_MSHPOLY_READ_ERROR_ELEMS_CEN_CLOSE_SUBSECTION:
		*msg = "Elements centroids subsection closing is not defined";
		break;
	case NB_MSHPOLY_READ_ERROR_ELEMS_NSIDES_OPEN_SUBSECTION:
		*msg = "Elements N-sides subsection opening is not defined";
		break;
	case NB_MSHPOLY_READ_ERROR_ELEMS_NSIDES_BAD_INPUT:
		*msg = "Bad N-sides definition in elements N-sides subsection";
		break;
	case NB_MSHPOLY_READ_ERROR_ELEMS_NSIDES_OUT_RANGE:
		*msg = "N-sides is out of range in elements N-sides subsection";
		break;
	case NB_MSHPOLY_READ_ERROR_ELEMS_NSIDES_CLOSE_SUBSECTION:
		*msg = "Elements N-sides subsection closing is not defined";
		break;
	case NB_MSHPOLY_READ_ERROR_ELEMS_ADJ_OPEN_SUBSECTION:
		*msg = "Elements adjacencies subsection opening is not defined";
		break;
	case NB_MSHPOLY_READ_ERROR_ELEMS_ADJ_BAD_INPUT:
		*msg = "Bad node ID definition in elements adj. subsection";
		break;
	case NB_MSHPOLY_READ_ERROR_ELEMS_ADJ_OUT_RANGE:
		*msg = "Node ID is out of range in elements adj. subsection";
		break;
	case NB_MSHPOLY_READ_ERROR_ELEMS_ADJ_CLOSE_SUBSECTION:
		*msg = "Elements adjacencies subsection closing is not defined";
		break;
	case NB_MSHPOLY_READ_ERROR_ELEMS_NGB_OPEN_SUBSECTION:
		*msg = "Elements neighbors subsection opening is not defined";
		break;
	case NB_MSHPOLY_READ_ERROR_ELEMS_NGB_BAD_INPUT:
		*msg = "Bad elem. ID definition in elem. neighbors subsection";
		break;
	case NB_MSHPOLY_READ_ERROR_ELEMS_NGB_OUT_RANGE:
		*msg = "Elem. ID is out of range in elem. neighbor subsection";
		break;
	case NB_MSHPOLY_READ_ERROR_ELEMS_NGB_CLOSE_SUBSECTION:
		*msg = "Elements neighbors subsection closing is not defined";
		break;
	case NB_MSHPOLY_READ_ERROR_ELEMS_CLOSE_SECTION:
		*msg = "Elements section closing is not defined";
		break;
	case NB_MSHPOLY_READ_ERROR_INVTX_OPEN_SECTION:
		*msg = "Input vertices section opening is not defined";
		break;
	case NB_MSHPOLY_READ_ERROR_INVTX_NOT_DEFINED:
		*msg = "Input vertices section exists but "	\
			NINVTX_VARNAME				\
			" = 0";
		break;
	case NB_MSHPOLY_READ_ERROR_INVTX_BAD_INPUT:
		*msg = "Bad node ID definition in input vertices section";
		break;
	case NB_MSHPOLY_READ_ERROR_INVTX_OUT_OF_RANGE:
		*msg = "Node ID is out of range in input vertices section";
		break;
	case NB_MSHPOLY_READ_ERROR_INVTX_CLOSE_SECTION:
		*msg = "Input vertices section closing is not defined";
		break;
	case NB_MSHPOLY_READ_ERROR_INSGM_OPEN_SECTION:
		*msg = "Input segments section opening is not defined";
		break;
	case NB_MSHPOLY_READ_ERROR_INSGM_NOT_DEFINED:
		*msg = "Input segment section exists but "	\
			NINSGM_VARNAME				\
			" = 0";
		break;
	case NB_MSHPOLY_READ_ERROR_INSGM_SIZE_OPEN_SUBSECTION:
		*msg = "Input segments size subsection opening is not defined";
		break;
	case NB_MSHPOLY_READ_ERROR_INSGM_SIZE_BAD_INPUT:
		*msg = "Bad size definition in input sgm size subsection";
		break;
	case NB_MSHPOLY_READ_ERROR_INSGM_SIZE_CLOSE_SUBSECTION:
		*msg = "Input segments size subsection closing is not defined";
		break;
	case NB_MSHPOLY_READ_ERROR_INSGM_IDS_OPEN_SUBSECTION:
		*msg = "Input segments IDs subsection opening is not defined";
		break;
	case NB_MSHPOLY_READ_ERROR_INSGM_IDS_OUT_OF_RANGE:
		*msg = "Node ID is out of range in input sgm IDs subsection";
		break;
	case NB_MSHPOLY_READ_ERROR_INSGM_IDS_BAD_OPEN_BRK:
		*msg = "List opening is not def. in input sgm IDs subsection";
		break;
	case NB_MSHPOLY_READ_ERROR_INSGM_IDS_BAD_INPUT:
		*msg = "Bad node ID definition in input sgm IDs subsection";
		break;
	case NB_MSHPOLY_READ_ERROR_INSGM_IDS_BAD_LINK:
		*msg = "Bad link definition in input sgm IDs subsection";
		break;
	case NB_MSHPOLY_READ_ERROR_INSGM_IDS_BAD_CLOSE_BRK:
		*msg = "List closing is not def. in input sgm IDs subsection";
		break;
	case NB_MSHPOLY_READ_ERROR_INSGM_IDS_CLOSE_SUBSECTION:
		*msg = "Input segments IDs subsection closing is not defined";
		break;
	case NB_MSHPOLY_READ_ERROR_INSGM_CLOSE_SECTION:
		*msg = "Input segments section closing is not defined";
		break;
	}
}
