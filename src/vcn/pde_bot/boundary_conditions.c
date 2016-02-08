#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>

#include "vcn/cfreader_cat.h"
#include "vcn/geometric_bot/utils2D.h"
#include "vcn/geometric_bot/mesh/elements2D/triangles.h"
#include "vcn/pde_bot/boundary_conditions.h"

vcn_bcond_t* vcn_fem_bcond_create(void)
{
	vcn_bcond_t* bconditions =
		(vcn_bcond_t*)calloc(1, sizeof(vcn_bcond_t));
	return bconditions;
}

inline vcn_bcond_t* vcn_fem_bcond_clone(const vcn_bcond_t*const bconditions)
{
	vcn_bcond_t *clone = calloc(1, sizeof(*clone));
	vcn_fem_bcond_copy(clone, bconditions);
	return clone;
}

vcn_bcond_t* vcn_fem_bcond_read(const char* filename)
{
	vcn_bcond_t* bconditions =
		(vcn_bcond_t*)calloc(1, sizeof(vcn_bcond_t));
  
	/* Initialize custom format to read file */
	vcn_cfreader_t* cfreader = vcn_cfreader_create(filename, "#");
	if (cfreader == NULL) {
		vcn_fem_bcond_destroy(bconditions);
		return NULL;
	}
	/* Read number of Dirichlet conditions upon vertices */
	if (vcn_cfreader_read_uint32_t(cfreader, 
				       &(bconditions->N_Dirichlet_on_vtx)) != 0) {
		vcn_cfreader_destroy(cfreader);
		vcn_fem_bcond_destroy(bconditions);
		return NULL;
	}
	/* Allocate data to store Dirichlet Cnd. upon vertices */
	bconditions->Dirichlet_on_vtx_idx = 
		malloc(bconditions->N_Dirichlet_on_vtx * sizeof(uint32_t));
	bconditions->Dirichlet_on_vtx_dof_mask =
		malloc(2 * bconditions->N_Dirichlet_on_vtx * sizeof(bool));
	bconditions->Dirichlet_on_vtx_val =
		calloc(2 * bconditions->N_Dirichlet_on_vtx, sizeof(double));

	/* Read number of Neuman conditions upon vertices */
	if (vcn_cfreader_read_uint32_t(cfreader,
				       &(bconditions->N_Neuman_on_vtx)) != 0) {
		vcn_cfreader_destroy(cfreader);
		vcn_fem_bcond_destroy(bconditions);
		return NULL;
	}
	/* Allocate data to store Neuman Cnd. upon vertices */
	bconditions->Neuman_on_vtx_idx = 
		(uint32_t*)malloc(bconditions->N_Neuman_on_vtx * sizeof(uint32_t));
	bconditions->Neuman_on_vtx_dof_mask =
		(bool*)malloc(2 * bconditions->N_Neuman_on_vtx * sizeof(bool));
	bconditions->Neuman_on_vtx_val =
		(double*)calloc(2 * bconditions->N_Neuman_on_vtx, sizeof(double));

	/* Read number of Dirichlet conditions upon segments */
	if (vcn_cfreader_read_uint32_t(cfreader,
				       &(bconditions->N_Dirichlet_on_sgm)) != 0) {
		vcn_cfreader_destroy(cfreader);
		vcn_fem_bcond_destroy(bconditions);
		return NULL;
	}
	/* Allocate data to store Dirichlet Cnd. upon segments */
	bconditions->Dirichlet_on_sgm_idx = 
		malloc(bconditions->N_Dirichlet_on_sgm * sizeof(uint32_t));
	bconditions->Dirichlet_on_sgm_dof_mask =
		malloc(2 * bconditions->N_Dirichlet_on_sgm * sizeof(bool));
	bconditions->Dirichlet_on_sgm_val =
		calloc(2 * bconditions->N_Dirichlet_on_sgm, sizeof(double));

	/* Read number of Neuman conditions upon segments */
	if (vcn_cfreader_read_uint32_t(cfreader, &(bconditions->N_Neuman_on_sgm)) != 0) {
		vcn_cfreader_destroy(cfreader);
		vcn_fem_bcond_destroy(bconditions);
		return NULL;
	}
	/* Allocate data to store Neuman Cnd. upon vertices */
	bconditions->Neuman_on_sgm_idx = 
		malloc(bconditions->N_Neuman_on_sgm * sizeof(uint32_t));
	bconditions->Neuman_on_sgm_dof_mask =
		malloc(2 * bconditions->N_Neuman_on_sgm * sizeof(bool));
	bconditions->Neuman_on_sgm_val =
		calloc(2 * bconditions->N_Neuman_on_sgm, sizeof(double));

	/* Read Dirichlet Cnd. upon vertices */
	for(uint32_t i = 0; i < bconditions->N_Dirichlet_on_vtx; i++){
		if (vcn_cfreader_read_uint32_t(cfreader, &(bconditions->Dirichlet_on_vtx_idx[i])) != 0) {
			vcn_cfreader_destroy(cfreader);
			vcn_fem_bcond_destroy(bconditions);
			return NULL;
		}
		uint32_t aux;
		if (vcn_cfreader_read_uint32_t(cfreader, &aux) != 0) {
			vcn_cfreader_destroy(cfreader);
			vcn_fem_bcond_destroy(bconditions);
			return NULL;
		}
		bconditions->Dirichlet_on_vtx_dof_mask[i * 2] = (aux == 1)?true:false;

		if (vcn_cfreader_read_uint32_t(cfreader, &aux) != 0) {
			vcn_cfreader_destroy(cfreader);
			vcn_fem_bcond_destroy(bconditions);
			return NULL;
		}
		bconditions->Dirichlet_on_vtx_dof_mask[i*2+1] = (aux == 1)?true:false;

		if (bconditions->Dirichlet_on_vtx_dof_mask[i * 2]) {
			if (vcn_cfreader_read_double(cfreader, &(bconditions->Dirichlet_on_vtx_val[i * 2])) != 0) {
				vcn_cfreader_destroy(cfreader);
				vcn_fem_bcond_destroy(bconditions);
				return NULL;
			}
		}

		if (bconditions->Dirichlet_on_vtx_dof_mask[i*2+1]) {
			if (vcn_cfreader_read_double(cfreader, &(bconditions->Dirichlet_on_vtx_val[i*2+1])) != 0) {
				vcn_cfreader_destroy(cfreader);
				vcn_fem_bcond_destroy(bconditions);
				return NULL;
			}
		}
	}

	/* Read Neuman Cnd. upon vertices */
	for (uint32_t i = 0; i < bconditions->N_Neuman_on_vtx; i++) {
		if (vcn_cfreader_read_uint32_t(cfreader, &(bconditions->Neuman_on_vtx_idx[i])) != 0) {
			vcn_cfreader_destroy(cfreader);
			vcn_fem_bcond_destroy(bconditions);
			return NULL;
		}
		uint32_t aux;
		if (vcn_cfreader_read_uint32_t(cfreader, &aux) != 0) {
			vcn_cfreader_destroy(cfreader);
			vcn_fem_bcond_destroy(bconditions);
			return NULL;
		}
		bconditions->Neuman_on_vtx_dof_mask[i * 2] = (aux == 1)?true:false;

		if (vcn_cfreader_read_uint32_t(cfreader, &aux) != 0) {
			vcn_cfreader_destroy(cfreader);
			vcn_fem_bcond_destroy(bconditions);
			return NULL;
		}
		bconditions->Neuman_on_vtx_dof_mask[i*2+1] = (aux == 1)?true:false;

		if (bconditions->Neuman_on_vtx_dof_mask[i * 2]) {
			if (vcn_cfreader_read_double(cfreader, &(bconditions->Neuman_on_vtx_val[i * 2])) != 0) {
				vcn_cfreader_destroy(cfreader);
				vcn_fem_bcond_destroy(bconditions);
				return NULL;
			}
		}

		if (bconditions->Neuman_on_vtx_dof_mask[i*2+1]) {
			if (vcn_cfreader_read_double(cfreader, &(bconditions->Neuman_on_vtx_val[i*2+1])) != 0) {
				vcn_cfreader_destroy(cfreader);
				vcn_fem_bcond_destroy(bconditions);
				return NULL;
			}
		}
	}

	/* Read Dirichlet Cnd. upon segments */
	for (uint32_t i = 0; i < bconditions->N_Dirichlet_on_sgm; i++) {
		if (vcn_cfreader_read_uint32_t(cfreader, &(bconditions->Dirichlet_on_sgm_idx[i])) != 0) {
			vcn_cfreader_destroy(cfreader);
			vcn_fem_bcond_destroy(bconditions);
			return NULL;
		}
		uint32_t aux;
		if (vcn_cfreader_read_uint32_t(cfreader, &aux) != 0) {
			vcn_cfreader_destroy(cfreader);
			vcn_fem_bcond_destroy(bconditions);
			return NULL;
		}
		bconditions->Dirichlet_on_sgm_dof_mask[i * 2] = (aux == 1)?true:false;

		if (vcn_cfreader_read_uint32_t(cfreader, &aux) != 0) {
			vcn_cfreader_destroy(cfreader);
			vcn_fem_bcond_destroy(bconditions);
			return NULL;
		}
		bconditions->Dirichlet_on_sgm_dof_mask[i*2+1] = (aux == 1)?true:false;

		if (bconditions->Dirichlet_on_sgm_dof_mask[i * 2]) {
			if(vcn_cfreader_read_double(cfreader, &(bconditions->Dirichlet_on_sgm_val[i * 2])) != 0){
				vcn_cfreader_destroy(cfreader);
				vcn_fem_bcond_destroy(bconditions);
				return NULL;
			}
		}

		if (bconditions->Dirichlet_on_sgm_dof_mask[i*2+1]) {
			if (vcn_cfreader_read_double(cfreader, &(bconditions->Dirichlet_on_sgm_val[i*2+1])) != 0){
				vcn_cfreader_destroy(cfreader);
				vcn_fem_bcond_destroy(bconditions);
				return NULL;
			}
		}
	}

	/* Read Neuman Cnd. upon vertices */
	for (uint32_t i = 0; i < bconditions->N_Neuman_on_sgm; i++) {
		if (vcn_cfreader_read_uint32_t(cfreader,
					       &(bconditions->Neuman_on_sgm_idx[i])) != 0) {
			vcn_cfreader_destroy(cfreader);
			vcn_fem_bcond_destroy(bconditions);
			return NULL;
		}
		uint32_t aux;
		if (vcn_cfreader_read_uint32_t(cfreader, &aux) != 0) {
			vcn_cfreader_destroy(cfreader);
			vcn_fem_bcond_destroy(bconditions);
			return NULL;
		}
		bconditions->Neuman_on_sgm_dof_mask[i * 2] = (aux == 1)?true:false;

		if (vcn_cfreader_read_uint32_t(cfreader, &aux) != 0) {
			vcn_cfreader_destroy(cfreader);
			vcn_fem_bcond_destroy(bconditions);
			return NULL;
		}
		bconditions->Neuman_on_sgm_dof_mask[i*2+1] = (aux == 1)?true:false;

		if(bconditions->Neuman_on_sgm_dof_mask[i * 2]){
			if(vcn_cfreader_read_double(cfreader, &(bconditions->Neuman_on_sgm_val[i * 2])) != 0){
				vcn_cfreader_destroy(cfreader);
				vcn_fem_bcond_destroy(bconditions);
				return NULL;
			}
		}

		if(bconditions->Neuman_on_sgm_dof_mask[i*2+1]){
			if(vcn_cfreader_read_double(cfreader, &(bconditions->Neuman_on_sgm_val[i*2+1])) != 0){
				vcn_cfreader_destroy(cfreader);
				vcn_fem_bcond_destroy(bconditions);
				return NULL;
			}
		}
	}

	/* Free memory */
	vcn_cfreader_destroy(cfreader);

	/* Return bconditions */
	return bconditions;
}

void vcn_fem_bcond_save(const vcn_bcond_t *const bconditions, 
			const char* filename){
	FILE* fp = fopen(filename, "w");
	if(fp == NULL) return;
  
	fprintf(fp, "# Boundary conditions description \n");
	fprintf(fp, "# Row <- [id x_mask y_mask (x_value) (y_value)]\n\n");

	fprintf(fp, "%i # Number of Dirichlet conditions upon vertices\n",
		bconditions->N_Dirichlet_on_vtx);
	fprintf(fp, "%i # Number of Neuman conditions upon vertices\n",
		bconditions->N_Neuman_on_vtx);
	fprintf(fp, "%i # Number of Dirichlet conditions upon segments\n",
		bconditions->N_Dirichlet_on_sgm);
	fprintf(fp, "%i # Number of Neuman conditions upon segments\n\n",
		bconditions->N_Neuman_on_sgm);

	fprintf(fp, "# Dirichlet conditions upon vertices  \n");
	for(uint32_t i = 0; i < bconditions->N_Dirichlet_on_vtx; i++){
		fprintf(fp, "%i %i %i ",
			bconditions->Dirichlet_on_vtx_idx[i],
			((bconditions->Dirichlet_on_vtx_dof_mask[i * 2])?1:0),
			((bconditions->Dirichlet_on_vtx_dof_mask[i*2+1])?1:0));
		if(bconditions->Dirichlet_on_vtx_dof_mask[i * 2])
			fprintf(fp, "%e ", bconditions->Dirichlet_on_vtx_val[i * 2]);
		if(bconditions->Dirichlet_on_vtx_dof_mask[i*2+1])
			fprintf(fp, "%e ", bconditions->Dirichlet_on_vtx_val[i*2+1]);
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n");

	fprintf(fp, "# Neuman conditions upon vertices  \n");
	for(uint32_t i = 0; i < bconditions->N_Neuman_on_vtx; i++){
		fprintf(fp, "%i %i %i ",
			bconditions->Neuman_on_vtx_idx[i],
			((bconditions->Neuman_on_vtx_dof_mask[i * 2])?1:0),
			((bconditions->Neuman_on_vtx_dof_mask[i*2+1])?1:0));
		if(bconditions->Neuman_on_vtx_dof_mask[i * 2])
			fprintf(fp, "%e ", bconditions->Neuman_on_vtx_val[i * 2]);
		if(bconditions->Neuman_on_vtx_dof_mask[i*2+1])
			fprintf(fp, "%e ", bconditions->Neuman_on_vtx_val[i*2+1]);
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n");

	fprintf(fp, "# Dirichlet conditions upon segments  \n");
	for(uint32_t i = 0; i < bconditions->N_Dirichlet_on_sgm; i++){
		fprintf(fp, "%i %i %i ",
			bconditions->Dirichlet_on_sgm_idx[i],
			((bconditions->Dirichlet_on_sgm_dof_mask[i * 2])?1:0),
			((bconditions->Dirichlet_on_sgm_dof_mask[i*2+1])?1:0));
		if(bconditions->Dirichlet_on_sgm_dof_mask[i * 2])
			fprintf(fp, "%e ", bconditions->Dirichlet_on_sgm_val[i * 2]);
		if(bconditions->Dirichlet_on_sgm_dof_mask[i*2+1])
			fprintf(fp, "%e ", bconditions->Dirichlet_on_sgm_val[i*2+1]);
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n");

	fprintf(fp, "# Neuman conditions upon segments  \n");
	for(uint32_t i = 0; i < bconditions->N_Neuman_on_sgm; i++){
		fprintf(fp, "%i %i %i ",
			bconditions->Neuman_on_sgm_idx[i],
			((bconditions->Neuman_on_sgm_dof_mask[i * 2])?1:0),
			((bconditions->Neuman_on_sgm_dof_mask[i*2+1])?1:0));
		if(bconditions->Neuman_on_sgm_dof_mask[i * 2])
			fprintf(fp, "%e ", bconditions->Neuman_on_sgm_val[i * 2]);
		if(bconditions->Neuman_on_sgm_dof_mask[i*2+1])
			fprintf(fp, "%e ", bconditions->Neuman_on_sgm_val[i*2+1]);
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n");  

	fclose(fp);
}

void vcn_fem_bcond_copy(vcn_bcond_t* cpy,
			const vcn_bcond_t *const src){
	cpy->N_dof = src->N_dof;
	cpy->N_Dirichlet_on_vtx = src->N_Dirichlet_on_vtx;
	cpy->N_Neuman_on_vtx = src->N_Neuman_on_vtx;
	cpy->N_Dirichlet_on_sgm = src->N_Dirichlet_on_sgm;
	cpy->N_Neuman_on_sgm = src->N_Neuman_on_sgm;

	if(cpy->N_Dirichlet_on_vtx > 0){
		cpy->Dirichlet_on_vtx_idx =
			(uint32_t*) malloc(cpy->N_Dirichlet_on_vtx * sizeof(uint32_t));
		cpy->Dirichlet_on_vtx_dof_mask =
			(bool*) malloc(cpy->N_dof * cpy->N_Dirichlet_on_vtx * sizeof(bool));
		cpy->Dirichlet_on_vtx_val =
			(double*) malloc(cpy->N_dof * cpy->N_Dirichlet_on_vtx * sizeof(double));

		memcpy(cpy->Dirichlet_on_vtx_idx, src->Dirichlet_on_vtx_idx,
		       cpy->N_Dirichlet_on_vtx * sizeof(uint32_t));
		memcpy(cpy->Dirichlet_on_vtx_dof_mask, src->Dirichlet_on_vtx_dof_mask,
		       cpy->N_dof * cpy->N_Dirichlet_on_vtx * sizeof(bool));
		memcpy(cpy->Dirichlet_on_vtx_val, src->Dirichlet_on_vtx_val,
		       cpy->N_dof * cpy->N_Dirichlet_on_vtx * sizeof(double));
	}

	if(cpy->N_Neuman_on_vtx > 0){
		cpy->Neuman_on_vtx_idx =
			(uint32_t*) malloc(cpy->N_Neuman_on_vtx * sizeof(uint32_t));
		cpy->Neuman_on_vtx_dof_mask =
			(bool*) malloc(cpy->N_dof * cpy->N_Neuman_on_vtx * sizeof(bool));
		cpy->Neuman_on_vtx_val =
			(double*) malloc(cpy->N_dof * cpy->N_Neuman_on_vtx * sizeof(double));

		memcpy(cpy->Neuman_on_vtx_idx, src->Neuman_on_vtx_idx,
		       cpy->N_Neuman_on_vtx * sizeof(uint32_t));
		memcpy(cpy->Neuman_on_vtx_dof_mask, src->Neuman_on_vtx_dof_mask,
		       cpy->N_dof * cpy->N_Neuman_on_vtx * sizeof(bool));
		memcpy(cpy->Neuman_on_vtx_val, src->Neuman_on_vtx_val,
		       cpy->N_dof * cpy->N_Neuman_on_vtx * sizeof(double));
	}

	if(cpy->N_Dirichlet_on_sgm > 0){
		cpy->Dirichlet_on_sgm_idx =
			(uint32_t*) malloc(cpy->N_Dirichlet_on_sgm * sizeof(uint32_t));
		cpy->Dirichlet_on_sgm_dof_mask =
			(bool*) malloc(cpy->N_dof * cpy->N_Dirichlet_on_sgm * sizeof(bool));
		cpy->Dirichlet_on_sgm_val =
			(double*) malloc(cpy->N_dof * cpy->N_Dirichlet_on_sgm * sizeof(double));

		memcpy(cpy->Dirichlet_on_sgm_idx, src->Dirichlet_on_sgm_idx,
		       cpy->N_Dirichlet_on_sgm * sizeof(uint32_t));
		memcpy(cpy->Dirichlet_on_sgm_dof_mask, src->Dirichlet_on_sgm_dof_mask,
		       cpy->N_dof * cpy->N_Dirichlet_on_sgm * sizeof(bool));
		memcpy(cpy->Dirichlet_on_sgm_val, src->Dirichlet_on_sgm_val,
		       cpy->N_dof * cpy->N_Dirichlet_on_sgm * sizeof(double));
	}

	if(cpy->N_Neuman_on_sgm > 0){
		cpy->Neuman_on_sgm_idx =
			(uint32_t*) malloc(cpy->N_Neuman_on_sgm * sizeof(uint32_t));
		cpy->Neuman_on_sgm_dof_mask =
			(bool*) malloc(cpy->N_dof * cpy->N_Neuman_on_sgm * sizeof(bool));
		cpy->Neuman_on_sgm_val =
			(double*) malloc(cpy->N_dof * cpy->N_Neuman_on_sgm * sizeof(double));

		memcpy(cpy->Neuman_on_sgm_idx, src->Neuman_on_sgm_idx,
		       cpy->N_Neuman_on_sgm * sizeof(uint32_t));
		memcpy(cpy->Neuman_on_sgm_dof_mask, src->Neuman_on_sgm_dof_mask,
		       cpy->N_dof * cpy->N_Neuman_on_sgm * sizeof(bool));
		memcpy(cpy->Neuman_on_sgm_val, src->Neuman_on_sgm_val,
		       cpy->N_dof * cpy->N_Neuman_on_sgm * sizeof(double));
	}
}

void vcn_fem_bcond_clear(vcn_bcond_t* bconditions){
	if(bconditions->N_Dirichlet_on_vtx > 0){
		free(bconditions->Dirichlet_on_vtx_idx);
		free(bconditions->Dirichlet_on_vtx_dof_mask);
		free(bconditions->Dirichlet_on_vtx_val);
	}
	if(bconditions->N_Neuman_on_vtx > 0){
		free(bconditions->Neuman_on_vtx_idx);
		free(bconditions->Neuman_on_vtx_dof_mask);
		free(bconditions->Neuman_on_vtx_val);
	}
	if(bconditions->N_Dirichlet_on_sgm > 0){
		free(bconditions->Dirichlet_on_sgm_idx);
		free(bconditions->Dirichlet_on_sgm_dof_mask);
		free(bconditions->Dirichlet_on_sgm_val);
	}
	if(bconditions->N_Neuman_on_sgm > 0){
		free(bconditions->Neuman_on_sgm_idx);
		free(bconditions->Neuman_on_sgm_dof_mask);
		free(bconditions->Neuman_on_sgm_val);
	}
	memset(bconditions, 0, sizeof(vcn_bcond_t));
}

void vcn_fem_bcond_destroy(vcn_bcond_t* bconditions){
	vcn_fem_bcond_clear(bconditions);
	free(bconditions);
}

void vcn_fem_bcond_printf(const vcn_bcond_t* const bconditions){
	printf("Dirichlet conditions on vertices: %i \n", bconditions->N_Dirichlet_on_vtx);
	for(uint32_t i = 0; i < bconditions->N_Dirichlet_on_vtx; i++){
		printf("  %i (%i %i) <- %e %e\n",
		       bconditions->Dirichlet_on_vtx_idx[i],
		       bconditions->Dirichlet_on_vtx_dof_mask[i * 2],
		       bconditions->Dirichlet_on_vtx_dof_mask[i*2+1],
		       bconditions->Dirichlet_on_vtx_val[i * 2],
		       bconditions->Dirichlet_on_vtx_val[i*2+1]);
	}
	printf("Dirichlet conditions on segments: %i \n", bconditions->N_Dirichlet_on_sgm);
	for(uint32_t i = 0; i < bconditions->N_Dirichlet_on_sgm; i++){
		printf("  %i (%i %i) <- %e %e\n",
		       bconditions->Dirichlet_on_sgm_idx[i],
		       bconditions->Dirichlet_on_sgm_dof_mask[i * 2],
		       bconditions->Dirichlet_on_sgm_dof_mask[i*2+1],
		       bconditions->Dirichlet_on_sgm_val[i * 2],
		       bconditions->Dirichlet_on_sgm_val[i*2+1]);
	}
	printf("Neumann conditions on vertices: %i \n", bconditions->N_Neuman_on_vtx);
	for(uint32_t i = 0; i < bconditions->N_Neuman_on_vtx; i++){
		printf("  %i (%i %i) <- %e %e\n",
		       bconditions->Neuman_on_vtx_idx[i],
		       bconditions->Neuman_on_vtx_dof_mask[i * 2],
		       bconditions->Neuman_on_vtx_dof_mask[i*2+1],
		       bconditions->Neuman_on_vtx_val[i * 2],
		       bconditions->Neuman_on_vtx_val[i*2+1]);
	}
	printf("Neuman conditions on segments: %i \n", bconditions->N_Neuman_on_sgm);
	for(uint32_t i = 0; i < bconditions->N_Neuman_on_sgm; i++){
		printf("  %i (%i %i) <- %e %e\n",
		       bconditions->Neuman_on_sgm_idx[i],
		       bconditions->Neuman_on_sgm_dof_mask[i * 2],
		       bconditions->Neuman_on_sgm_dof_mask[i*2+1],
		       bconditions->Neuman_on_sgm_val[i * 2],
		       bconditions->Neuman_on_sgm_val[i*2+1]);
	}
}

vcn_bcond_t* vcn_fem_bcond_create_from_model_to_mesh
(const vcn_msh3trg_t *const msh3trg,
 const vcn_bcond_t *const bcond_on_model)
/* Mapping of the boundary conditions on the model to the boundary 
 * conditions on the mesh. The conditions are placed on the nodes.
 */
{
	vcn_bcond_t* bmeshcond = calloc(1, sizeof(*bmeshcond));

	bmeshcond->N_dof = bcond_on_model->N_dof;
	/***** Get Dirichlet conditions *****/
	bmeshcond->N_Dirichlet_on_vtx = bcond_on_model->N_Dirichlet_on_vtx;
	/* Set Dirichlet conditions on vertices forming the segments with 
	 * the original conditions.
	 */
	for (uint32_t i = 0; i < bcond_on_model->N_Dirichlet_on_sgm; i++) {
		uint32_t geom_idx = bcond_on_model->Dirichlet_on_sgm_idx[i];
		if (0 == msh3trg->N_subsgm_x_inputsgm[geom_idx])
			continue;
		bmeshcond->N_Dirichlet_on_vtx += 
			msh3trg->N_subsgm_x_inputsgm[geom_idx] + 1;
	}
	/* Allocate Dirichlet conditions */
	if(bmeshcond->N_Dirichlet_on_vtx > 0){
		bmeshcond->Dirichlet_on_vtx_idx =
			calloc(bmeshcond->N_Dirichlet_on_vtx,
					  sizeof(uint32_t));
		bmeshcond->Dirichlet_on_vtx_dof_mask =
			malloc(bmeshcond->N_dof * bmeshcond->N_Dirichlet_on_vtx *
				      sizeof(bool));
		bmeshcond->Dirichlet_on_vtx_val =
			malloc(bmeshcond->N_dof * bmeshcond->N_Dirichlet_on_vtx *
					sizeof(double));

		/* Set Dirichlet conditions on vertices */
		for(uint32_t i = 0; i < bcond_on_model->N_Dirichlet_on_vtx; i++){
			uint32_t geom_idx = bcond_on_model->Dirichlet_on_vtx_idx[i];
			bmeshcond->Dirichlet_on_vtx_idx[i] =
				msh3trg->input_vertices[geom_idx];
		}

		memcpy(bmeshcond->Dirichlet_on_vtx_dof_mask,
		       bcond_on_model->Dirichlet_on_vtx_dof_mask, sizeof(bool) *
		       bcond_on_model->N_Dirichlet_on_vtx * bcond_on_model->N_dof);
		memcpy(bmeshcond->Dirichlet_on_vtx_val,
		       bcond_on_model->Dirichlet_on_vtx_val, sizeof(double) *
		       bcond_on_model->N_Dirichlet_on_vtx * bcond_on_model->N_dof);

		/* Set Dirichlet conditions on segments */
		uint32_t idx = bcond_on_model->N_Dirichlet_on_vtx;
		for(uint32_t i=0; i < bcond_on_model->N_Dirichlet_on_sgm; i++){
			uint32_t geom_idx = bcond_on_model->Dirichlet_on_sgm_idx[i];
			if(msh3trg->N_subsgm_x_inputsgm[geom_idx] == 0) continue;
			for(uint32_t j = 0; j < msh3trg->N_subsgm_x_inputsgm[geom_idx] + 1; j++){
				bmeshcond->Dirichlet_on_vtx_idx[idx] = 
					msh3trg->meshvtx_x_inputsgm[geom_idx][j];

				for(uint32_t k = 0; k < bcond_on_model->N_dof; k++){
					bmeshcond->Dirichlet_on_vtx_dof_mask[idx * bcond_on_model->N_dof + k] = 
						bcond_on_model->Dirichlet_on_sgm_dof_mask[i * bcond_on_model->N_dof + k];
				}

				for(uint32_t k = 0; k < bcond_on_model->N_dof; k++){	
					bmeshcond->Dirichlet_on_vtx_val[idx * bcond_on_model->N_dof + k] = 
						bcond_on_model->Dirichlet_on_sgm_val[i * bcond_on_model->N_dof + k];
				}

				idx ++;
			}
		}
	}

	/***** Get Neuman conditions *****/
	bmeshcond->N_Neuman_on_vtx = bcond_on_model->N_Neuman_on_vtx;
	/* Set Neuman conditions on vertices forming the segments with 
	 * the original conditions.
	 */
	for(uint32_t i=0; i < bcond_on_model->N_Neuman_on_sgm; i++){
		uint32_t geom_idx = bcond_on_model->Neuman_on_sgm_idx[i];
		if(msh3trg->N_subsgm_x_inputsgm[geom_idx] == 0) continue;
		bmeshcond->N_Neuman_on_vtx += 
			msh3trg->N_subsgm_x_inputsgm[geom_idx] + 1;
	}
	/* Allocate Neuman conditions */
	if(bmeshcond->N_Neuman_on_vtx > 0){
		bmeshcond->Neuman_on_vtx_idx =
			(uint32_t*)calloc(bmeshcond->N_Neuman_on_vtx, sizeof(uint32_t));
		bmeshcond->Neuman_on_vtx_dof_mask =
			(bool*)malloc(bmeshcond->N_dof * bmeshcond->N_Neuman_on_vtx *
				      sizeof(bool));
		bmeshcond->Neuman_on_vtx_val =
			(double*)malloc(bmeshcond->N_dof * bmeshcond->N_Neuman_on_vtx *
					sizeof(double));
    
		/* Set Neuman conditions on vertices */
		for(uint32_t i=0; i < bcond_on_model->N_Neuman_on_vtx; i++){
			uint32_t geom_idx = bcond_on_model->Neuman_on_vtx_idx[i];
			bmeshcond->Neuman_on_vtx_idx[i] =
				msh3trg->input_vertices[geom_idx];
		}

		memcpy(bmeshcond->Neuman_on_vtx_dof_mask,
		       bcond_on_model->Neuman_on_vtx_dof_mask, sizeof(bool) *
		       bcond_on_model->N_Neuman_on_vtx * bcond_on_model->N_dof);
		memcpy(bmeshcond->Neuman_on_vtx_val,
		       bcond_on_model->Neuman_on_vtx_val, sizeof(double) *
		       bcond_on_model->N_Neuman_on_vtx * bcond_on_model->N_dof);

		/* Set Neuman conditions on segments */
		uint32_t idx = bcond_on_model->N_Neuman_on_vtx;
		for(uint32_t i=0; i < bcond_on_model->N_Neuman_on_sgm; i++){
			uint32_t geom_idx = bcond_on_model->Neuman_on_sgm_idx[i];
			if(msh3trg->N_subsgm_x_inputsgm[geom_idx] == 0) continue;
      
			uint32_t s1 = msh3trg->meshvtx_x_inputsgm[geom_idx][0];
			uint32_t s2 = msh3trg->meshvtx_x_inputsgm[geom_idx]
				[msh3trg->N_subsgm_x_inputsgm[geom_idx]];
			double sgm_length = vcn_utils2D_get_dist(&(msh3trg->vertices[s1 * 2]),
								 &(msh3trg->vertices[s2 * 2]));

			for(uint32_t j = 0; j < msh3trg->N_subsgm_x_inputsgm[geom_idx] + 1; j++){
				uint32_t sj = msh3trg->meshvtx_x_inputsgm[geom_idx][j];
				bmeshcond->Neuman_on_vtx_idx[idx] = sj;

				for(uint32_t k = 0; k < bcond_on_model->N_dof; k++){
					bmeshcond->Neuman_on_vtx_dof_mask[idx * bcond_on_model->N_dof + k] = 
						bcond_on_model->Neuman_on_sgm_dof_mask[i * bcond_on_model->N_dof + k];
				}

				/* Compute proportional weight to assign condition */
				double w_left = 0.0;
				if(j > 0){
					uint32_t sj_prev = msh3trg->meshvtx_x_inputsgm[geom_idx][j-1];
					double dist =
						vcn_utils2D_get_dist(&(msh3trg->vertices[sj_prev * 2]),
								     &(msh3trg->vertices[sj * 2]));
					w_left = 0.5 * (dist);
				}

				double w_right = 0.0;
				if(j < msh3trg->N_subsgm_x_inputsgm[geom_idx]){
					uint32_t sj_next = msh3trg->meshvtx_x_inputsgm[geom_idx][j+1];
					double dist =
						vcn_utils2D_get_dist(&(msh3trg->vertices[sj_next * 2]),
								     &(msh3trg->vertices[sj * 2]));
					w_right = 0.5 * (dist);
				}
	
				double weight = (w_left + w_right) / sgm_length;

				for(uint32_t k = 0; k < bcond_on_model->N_dof; k++){	
					bmeshcond->Neuman_on_vtx_val[idx * bcond_on_model->N_dof + k] = 
						bcond_on_model->Neuman_on_sgm_val[i * bcond_on_model->N_dof + k] * weight;
				}

				idx ++;
			}
		}
	}

	return bmeshcond;
}
