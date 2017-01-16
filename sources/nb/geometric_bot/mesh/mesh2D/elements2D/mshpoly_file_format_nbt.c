#include <stdio.h>
#include <stdint.h>

#include "nb/cfreader_bot.h"
#include "nb/geometric_bot.h"

#include "mshpoly_struct.h"
#include "mshpoly_file_format_nbt.h"

static void write_nodes(FILE *fp, const void *msh_ptr);
static void write_edges(FILE *fp, const void *msh_ptr);
static void write_elems(FILE *fp, const void *msh_ptr);
static void write_invtx(FILE *fp, const void *msh_ptr);
static void write_insgm(FILE *fp, const void *msh_ptr);

void nb_mshpoly_write_data_nbt(FILE *fp, const void *msh)
{
	write_nodes(fp, msh);
	write_edges(fp, msh);
	write_elems(fp, msh);
	write_invtx(fp, msh);
	write_insgm(fp, msh);
}

static void write_nodes(FILE *fp, const void *msh_ptr)
{
	const nb_mshpoly_t *msh = msh_ptr;

	fprintf(fp, "N_nodes = %i\n", msh->N_nod);
	for (uint32_t i = 0; i < msh->N_nod; i++)
		fprintf(fp, "  %e  %e\n", msh->nod[i * 2], msh->nod[i*2+1]);
	fprintf(fp, "\n");
}

static void write_edges(FILE *fp, const void *msh_ptr)
{
	const nb_mshpoly_t *msh = msh_ptr;

	fprintf(fp, "N_edges = %i\n", msh->N_edg);
	for (uint32_t i = 0; i < msh->N_edg; i++)
		fprintf(fp, "  %i %i\n",
			msh->edg[i * 2] + 1,
			msh->edg[i*2+1] + 1);
	fprintf(fp, "\n");
}

static void write_elems(FILE *fp, const void *msh_ptr)
{
	const nb_mshpoly_t *msh = msh_ptr;

	fprintf(fp, "N_elems = %i\n", msh->N_elems);
	fprintf(fp, "  Centroids\n");
	for (uint32_t i = 0; i < msh->N_elems; i++)
		fprintf(fp, "    %e %e\n", msh->cen[i * 2], msh->cen[i*2+1]);
	fprintf(fp, "\n");

	fprintf(fp, "  Adjacencies\n");
	for (uint32_t i = 0; i < msh->N_elems; i++) {
		fprintf(fp, "    %i ", msh->N_adj[i]);
		for (uint32_t j = 0; j < msh->N_adj[i]; j++)
			fprintf(fp, "%i ", msh->adj[i][j] + 1);
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n");

	fprintf(fp, "  Neighbors\n");
	fprintf(fp, "  # The null neighbors at the boundary have an ID\n");
	fprintf(fp, "  # greater than the number of elements.\n");
	for (uint32_t i = 0; i < msh->N_elems; i++) {
		fprintf(fp, "    %i ", msh->N_adj[i]);
		for (uint32_t j = 0; j < msh->N_adj[i]; j++)
			fprintf(fp, "%i ", msh->ngb[i][j] + 1);
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n");
}

static void write_invtx(FILE *fp, const void *msh_ptr)
{
	const nb_mshpoly_t *msh = msh_ptr;

	fprintf(fp, "N_input_vtx = %i\n", msh->N_vtx);
	for (uint32_t i = 0; i < msh->N_vtx; i++)
		fprintf(fp, " %i ->  %i\n", i + 1, msh->vtx[i] + 1);
	fprintf(fp, "\n");
}

static void write_insgm(FILE *fp, const void *msh_ptr)
{
	const nb_mshpoly_t *msh = msh_ptr;

	fprintf(fp, "N_input_sgm = %i\n", msh->N_sgm);
	for (uint32_t i = 0; i < msh->N_sgm; i++) {
		fprintf(fp, "    %i", msh->N_nod_x_sgm[i]);
		for (uint32_t j = 0; j < msh->N_nod_x_sgm[i]; j++)
			fprintf(fp, "%i ", msh->nod_x_sgm[i][j] + 1);
	}
	fprintf(fp, "\n");
}

int nb_mshpoly_read_data_nbt(nb_cfreader_t *cfr, void *msh)
{
	return 1;
}
