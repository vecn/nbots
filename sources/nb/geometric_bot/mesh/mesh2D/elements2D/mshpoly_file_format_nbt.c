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
static int read_elems_ngb(nb_cfreader_t *cfr, nb_mshpoly_t *msh);
static int read_invtx(nb_cfreader_t *cfr, nb_mshpoly_t *msh,
		      bool *exist_section);
static int read_insgm(nb_cfreader_t *cfr, nb_mshpoly_t *msh,
		      bool *exist_section);
static int read_insgm_size(nb_cfreader_t *cfr, nb_mshpoly_t *msh);
static int read_insgm_idvtx(nb_cfreader_t *cfr, nb_mshpoly_t *msh);

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

	write_elems_centroids(fp, msh);
	fprintf(fp, "\n");

	write_elems_N_sides(fp, msh);
	fprintf(fp, "\n");

	write_elems_adj(fp, msh);
	fprintf(fp, "\n");

	write_elems_ngb(fp, msh);

	fprintf(fp, "End Elems\n");
}

static void write_elems_centroids(FILE *fp, const nb_mshpoly_t *msh)
{
	fprintf(fp, "  Centroids\n");
	for (uint32_t i = 0; i < msh->N_elems; i++)
		fprintf(fp, "    %e %e\n", msh->cen[i * 2], msh->cen[i*2+1]);
	fprintf(fp, "  End Centroids\n");
}

static void write_elems_N_sides(FILE *fp, const nb_mshpoly_t *msh)
{
	fprintf(fp, "  N_sides\n");
	for (uint32_t i = 0; i < msh->N_elems; i++) {
		if (0 == i % 10)
			fprintf(fp, "    ");
		fprintf(fp, "%i  ", msh->N_adj[i]);
		
		if (0 == (i+1) % 10 ||  (i+1) == msh->N_elems)
			fprintf(fp, "\n");
	}
	fprintf(fp, "  End N_sides\n");
}

static void write_elems_adj(FILE *fp, const nb_mshpoly_t *msh)
{
	fprintf(fp, "  Adjacencies\n");
	for (uint32_t i = 0; i < msh->N_elems; i++) {
		fprintf(fp, "    ");
		for (uint32_t j = 0; j < msh->N_adj[i]; j++)
			fprintf(fp, "%i ", msh->adj[i][j] + 1);
		fprintf(fp, "\n");
	}
	fprintf(fp, "  End Adjacencies\n");
}

static void write_elems_ngb(FILE *fp, const nb_mshpoly_t *msh)
{
	fprintf(fp, "  Neighbors\n");
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
	fprintf(fp, "  End Neighbors\n");
}

static void write_invtx(FILE *fp, const nb_mshpoly_t *msh)
{
	int N_null = 0;
	for (uint32_t i = 0; i < msh->N_vtx; i++) {
		if (msh->vtx[i] < msh->N_vtx)
			N_null += 1;
	}
	if (N_null < msh->N_vtx) {
		fprintf(fp, "Input vertices\n");
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
		fprintf(fp, "End input vertices\n");
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
		fprintf(fp, "Input segments\n");
		if (N_null > 0)
			fprintf(fp, "  # Size = 0 indicates that mesh" \
				"segment does not exist.\n");
		write_insgm_size(fp, msh);
		fprintf(fp, "\n");
		write_insgm_idvtx(fp, msh);
		fprintf(fp, "End input segments\n");
	}
}

static void write_insgm_size(FILE *fp, const nb_mshpoly_t *msh)
{
	fprintf(fp, "  N_size\n");
	for (uint32_t i = 0; i < msh->N_sgm; i++) {
		if (0 == i % 10)
			fprintf(fp, "    ");
		fprintf(fp, "%i ", msh->N_nod_x_sgm[i]);
		
		if (0 == (i+1) % 10 ||  (i+1) == msh->N_elems)
			fprintf(fp, "\n");
	}
	fprintf(fp, "  End N_size\n\n");
}

static void write_insgm_idvtx(FILE *fp, const nb_mshpoly_t *msh)
{
	fprintf(fp, "  ID_vtx\n");
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
	fprintf(fp, "  End ID_vtx\n");
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

	bool exist_section;
	status = read_invtx(cfr, msh, &exist_section);
	if (0 != status)
		goto EXIT;

	if (!exist_section && msh->N_vtx > 0) {
		status = 1;
		goto EXIT;
	}

	status = read_insgm(cfr, msh, &exist_section);
	if (0 != status)
		goto EXIT;

	if (!exist_section && msh->N_sgm > 0)
		status = 1;
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
	int status = nb_cfreader_read_var_uint(cfr, "N_nodes", &N);
	if (0 != status)
		goto EXIT;
	msh->N_nod = N;

	status = nb_cfreader_read_var_uint(cfr, "N_edges", &N);
	if (0 != status)
		goto EXIT;
	msh->N_edg = N;

	status = nb_cfreader_read_var_uint(cfr, "N_elems", &N);
	if (0 != status)
		goto EXIT;
	msh->N_elems = N;

	status = nb_cfreader_read_var_uint(cfr, "N_input_vtx", &N);
	if (0 != status)
		goto EXIT;
	msh->N_vtx = N;

	status = nb_cfreader_read_var_uint(cfr, "N_input_sgm", &N);
	if (0 != status)
		goto EXIT;
	msh->N_sgm = N;

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
		if (0 == val || val > msh->N_nod) {
			status = 1;
			goto EXIT;
		}
		msh->edg[i * 2] = val - 1;

		status = nb_cfreader_read_uint(cfr, &val);
		if (0 != status)
			goto EXIT;
		if (0 == val || val > msh->N_nod) {
			status = 1;
			goto EXIT;
		}
		msh->edg[i*2+1] = val - 1;
	}

	status = nb_cfreader_check_line(cfr, "End Edges");
EXIT:
	return status;
}

static int read_elems(nb_cfreader_t *cfr, nb_mshpoly_t *msh)
{
	int status = nb_cfreader_check_line(cfr, "Elems");
	if (0 != status)
		goto EXIT;
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

	status = nb_cfreader_check_line(cfr, "End Elems");
EXIT:
	return status;
}

static int read_elems_centroids(nb_cfreader_t *cfr, nb_mshpoly_t *msh)
{
	int status = nb_cfreader_check_line(cfr, "Centroids");
	if (0 != status)
		goto EXIT;

	for (uint32_t i = 0; i < msh->N_elems; i++) {
		double val;
		status = nb_cfreader_read_double(cfr, &val);
		if (0 != status)
			goto EXIT;
		msh->cen[i * 2] = val;

		status = nb_cfreader_read_double(cfr, &val);
		if (0 != status)
			goto EXIT;
		msh->cen[i*2+1] = val;
	}
	status = nb_cfreader_check_line(cfr, "End Centroids");
EXIT:
	return status;
}

static int read_elems_N_sides(nb_cfreader_t *cfr, nb_mshpoly_t *msh)
{
	int status = nb_cfreader_check_line(cfr, "N_sides");
	if (0 != status)
		goto EXIT;

	for (uint32_t i = 0; i < msh->N_elems; i++) {
		uint32_t val;
		status = nb_cfreader_read_uint(cfr, &val);
		if (0 != status)
			goto EXIT;
		if (0 == val) {
			status = 1;
			goto EXIT;
		}
		msh->N_adj[i] = val;
	}

	status = nb_cfreader_check_line(cfr, "End N_sides");
EXIT:
	return status;
}

static int read_elems_adj(nb_cfreader_t *cfr, nb_mshpoly_t *msh)
{
	int status = nb_cfreader_check_line(cfr, "Adjacencies");
	if (0 != status)
		goto EXIT;

	for (uint32_t i = 0; i < msh->N_elems; i++) {
		for (int j = 0; j < msh->N_adj[i]; j++) {
			uint32_t val;
			status = nb_cfreader_read_uint(cfr, &val);
			if (0 != status)
				goto EXIT;
			if (0 == val || val > msh->N_nod) {
				status = 1;
				goto EXIT;
			}
			msh->adj[i][j] = val - 1;
		}
	}

	status = nb_cfreader_check_line(cfr, "End Adjacencies");
EXIT:
	return status;
}

static int read_elems_ngb(nb_cfreader_t *cfr, nb_mshpoly_t *msh)
{
	int status = nb_cfreader_check_line(cfr, "Neighbors");
	if (0 != status)
		goto EXIT;

	for (uint32_t i = 0; i < msh->N_elems; i++) {
		for (int j = 0; j < msh->N_adj[i]; j++) {
			uint32_t val;
			status = nb_cfreader_read_uint(cfr, &val);
			if (0 != status)
				goto EXIT;
			if (0 == val) {
				msh->ngb[i][j] = msh->N_elems;
			} else if (val < msh->N_elems) {
				msh->ngb[i][j] = val - 1;
			} else {
				status = 1;
				goto EXIT;
			}
		}
	}

	status = nb_cfreader_check_line(cfr, "End Neighbors");
EXIT:
	return status;
}

static int read_invtx(nb_cfreader_t *cfr, nb_mshpoly_t *msh,
		      bool *exist_section)
{
	int status = nb_cfreader_check_line(cfr, "Input vertices");
	if (0 != status) {
		status = 0;
		*exist_section = false;
		goto EXIT;
	}
	*exist_section = true;

	for (uint32_t i = 0; i < msh->N_vtx; i++) {
		char val1[15];
		char val2[15];
		status = nb_cfreader_read_tuple(cfr, val1, val2);
		if (0 != status)
			goto EXIT;
		
		status = 1;
		char *pch;
		unsigned int N = strtoul(val1, &pch, 10);
		if (pch == val1)
			goto EXIT;

		if (0 == N || N > msh->N_vtx)
			goto EXIT;
		uint32_t id = N-1;

		N = strtoul(val2, &pch, 10);
		if (pch == val1)
			goto EXIT;

		if (N > 0 && N <= msh->N_nod)
			msh->vtx[id] = N - 1;
		else if (0 == N)
			msh->vtx[id] = msh->N_nod;
		else
			goto EXIT;
	}

	status = nb_cfreader_check_line(cfr, "End input vertices");
EXIT:
	return status;
}

static int read_insgm(nb_cfreader_t *cfr, nb_mshpoly_t *msh,
		      bool *exist_section)
{
	int status = nb_cfreader_check_line(cfr, "Input segments");
	if (0 != status) {
		status = 0;
		*exist_section = false;
		goto EXIT;
	}
	*exist_section = true;

	status = read_insgm_size(cfr, msh);
	if (0 != status)
		goto EXIT;

	nb_mshpoly_set_mem_of_nod_x_sgm(msh);

	status = read_insgm_idvtx(cfr, msh);
	if (0 != status)
		goto EXIT;

	status = nb_cfreader_check_line(cfr, "End input segments");
EXIT:
	return status;
}

static int read_insgm_size(nb_cfreader_t *cfr, nb_mshpoly_t *msh)
{
	int status = nb_cfreader_check_line(cfr, "N_size");
	if (0 != status)
		goto EXIT;

	for (uint32_t i = 0; i < msh->N_sgm; i++) {
		uint32_t val;
		status = nb_cfreader_read_uint(cfr, &val);
		if (0 != status)
			goto EXIT;

		msh->N_nod_x_sgm[i] = val;
	}

	status = nb_cfreader_check_line(cfr, "End N_size");
EXIT:
	return status;

}

static int read_insgm_idvtx(nb_cfreader_t *cfr, nb_mshpoly_t *msh)
{
	int status = nb_cfreader_check_line(cfr, "ID_vtx");
	if (0 != status)
		goto EXIT;

	for (uint32_t i = 0; i < msh->N_sgm; i++) {
		int status = nb_cfreader_check_token(cfr, "[");
		if (0 != status)
			goto EXIT;
		for (uint32_t j = 0; j < msh->N_nod_x_sgm[i]; j++) {
			if (j > 0) {
				status = nb_cfreader_check_token(cfr, "->");
				if (0 != status)
					goto EXIT;
			}
			uint32_t val;
			status = nb_cfreader_read_uint(cfr, &val);
			if (0 != status)
				goto EXIT;

			if (0 == val || val > msh->N_nod) {
				status = 1;
				goto EXIT;
			}

			msh->nod_x_sgm[i][j] = val - 1;
		}
		status = nb_cfreader_check_token(cfr, "]");
		if (0 != status)
			goto EXIT;
	}

	status = nb_cfreader_check_line(cfr, "End ID_vtx");
EXIT:
	return status;
}
