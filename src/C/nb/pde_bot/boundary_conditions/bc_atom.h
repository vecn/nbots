#ifndef __NB_PDE_BOT_BOUNDARY_CONDITIONS_BC_ATOM_H__
#define __NB_PDE_BOT_BOUNDARY_CONDITIONS_BC_ATOM_H__

typedef struct {
	uint32_t id;
	bool *mask;
	double *val;
} bc_atom_t;

bc_atom_t* bc_atom_create(uint8_t N_dof);
bc_atom_t* bc_atom_clone(const void *const bc_ptr, uint8_t N_dof);
void bc_atom_destroy(void *bc_ptr);

#endif
