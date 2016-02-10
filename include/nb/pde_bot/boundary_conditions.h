#ifndef __NB_PDE_BOT_BOUNDARY_CONDITIONS_H__
#define __NB_PDE_BOT_BOUNDARY_CONDITIONS_H__

#include "nb/geometric_bot/mesh/elements2D/triangles.h"

typedef struct {
	/* TEMPORAL: Change Neuman for Neumann */
	/* Boundary conditions */
	uint8_t N_dof; /* Degrees of freedom */
	/* Dirichlet conditions on vertices */
	uint32_t N_Dirichlet_on_vtx;
	uint32_t *Dirichlet_on_vtx_idx;
	bool *Dirichlet_on_vtx_dof_mask;
	double *Dirichlet_on_vtx_val;
	/* Neuman conditions on vertices */
	uint32_t N_Neuman_on_vtx;
	uint32_t *Neuman_on_vtx_idx;
	bool *Neuman_on_vtx_dof_mask;
	double *Neuman_on_vtx_val;
	/* Dirichlet conditions on segments */
	uint32_t N_Dirichlet_on_sgm;
	uint32_t *Dirichlet_on_sgm_idx;
	bool*Dirichlet_on_sgm_dof_mask;
	double *Dirichlet_on_sgm_val;
	/* Neuman conditions on segments */
	uint32_t N_Neuman_on_sgm;
	uint32_t *Neuman_on_sgm_idx;
	bool *Neuman_on_sgm_dof_mask;
	double *Neuman_on_sgm_val;
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

vcn_bcond_t* vcn_fem_bcond_create_from_model_to_mesh
			(const vcn_msh3trg_t *const msh3trg,
			 const vcn_bcond_t *const bconditions);

#endif
