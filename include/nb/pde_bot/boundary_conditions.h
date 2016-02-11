#ifndef __NB_PDE_BOT_BOUNDARY_CONDITIONS_H__
#define __NB_PDE_BOT_BOUNDARY_CONDITIONS_H__

#include "nb/geometric_bot/mesh/elements2D/triangles.h"

typedef struct {
	uint32_t N;
	uint32_t *id;
	bool *dof_mask;
	double *val;
} bc_data_t;

typedef struct {
	/* TEMPORAL: Change Neuman for Neumann */
	/* Boundary conditions */
	uint8_t N_dof; /* Degrees of freedom */
	bc_data_t dirichlet_on_vtx;
	bc_data_t neumann_on_vtx;
	bc_data_t dirichlet_on_sgm;
	bc_data_t neumann_on_sgm;
} vcn_bcond_t;

vcn_bcond_t* vcn_fem_bcond_create(void);
vcn_bcond_t* vcn_fem_bcond_clone(const vcn_bcond_t *const bconditions);
vcn_bcond_t* vcn_fem_bcond_read(const char* filename);
void vcn_fem_bcond_save(const vcn_bcond_t *const bconditions, 
			const char* filename);
void vcn_fem_bcond_copy(vcn_bcond_t* cpy,
			const vcn_bcond_t *const src);
void vcn_fem_bcond_clear(vcn_bcond_t* bconditions);
void vcn_fem_bcond_destroy(vcn_bcond_t* bconditions);
void vcn_fem_bcond_printf(const vcn_bcond_t* const bconditions);

#endif
