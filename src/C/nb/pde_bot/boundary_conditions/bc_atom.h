#ifndef __NB_PDE_BOT_BOUNDARY_CONDITIONS_BC_ATOM_H__
#define __NB_PDE_BOT_BOUNDARY_CONDITIONS_BC_ATOM_H__

#include <stdint.h>

typedef struct {
	uint32_t id;
	bool *mask;
	double *val;
} bc_atom_t;

uint16_t bc_atom_get_memsize(uint8_t N_dof);
bc_atom_t* bc_atom_create(uint8_t N_dof);
bc_atom_t* bc_atom_clone(const void *const bc_ptr, uint8_t N_dof);
void bc_atom_destroy(void *bc_ptr);
void bc_atom_set_data(bc_atom_t *bc,uint32_t elem_id,
		      const bool dof_mask[], const double value[]);

#endif
