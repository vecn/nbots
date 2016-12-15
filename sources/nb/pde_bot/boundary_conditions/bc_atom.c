#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "nb/memory_bot.h"
#include "bc_atom.h"

static void init_data_ptr(void *bc_ptr, uint8_t N_dof);
static void* allocate_bc_atom(uint8_t N_dof);

uint16_t bc_atom_get_memsize(uint8_t N_dof)
{
	uint16_t size = sizeof(bc_atom_t);
	uint16_t size_mask = sizeof(bool);
	uint16_t size_val = sizeof(double);
	return size + N_dof * (size_mask + size_val);
}

inline void bc_atom_init(void *bc_ptr, uint8_t N_dof)
{
	memset(bc_ptr, 0, bc_atom_get_memsize(N_dof));
	init_data_ptr(bc_ptr, N_dof);
}

static void init_data_ptr(void *bc_ptr, uint8_t N_dof)
{
	bc_atom_t *bc = bc_ptr;
	uint16_t size = sizeof(bc_atom_t);
	uint16_t size_mask = sizeof(*(bc->mask));
	char *memblock = bc_ptr;
	bc->mask = (void*)(memblock + size);
	bc->val = (void*)(memblock + size + N_dof * size_mask);
	bc->fval = NULL;
}

void bc_atom_copy(void *bc_ptr, const void *src_bc_ptr, uint8_t N_dof)
{
	init_data_ptr(bc_ptr, N_dof);
	bc_atom_t *bc = bc_ptr;
	const bc_atom_t *src_bc = src_bc_ptr;

	bc->id = src_bc->id;
	memcpy(bc->mask, src_bc->mask, N_dof * sizeof(*(bc->mask)));
	memcpy(bc->val, src_bc->val, N_dof * sizeof(*(bc->val)));
	bc->fval = src_bc->fval;
}

void bc_atom_finish(void *bc_ptr)
{
	; /* Null statement */
}

void* bc_atom_create(uint8_t N_dof)
{
	void *bc = allocate_bc_atom(N_dof);
	bc_atom_init(bc, N_dof);
	return bc;
}

static inline void* allocate_bc_atom(uint8_t N_dof)
{
	uint16_t size = bc_atom_get_memsize(N_dof);
	return nb_allocate_mem(size);
}

void* bc_atom_clone(const void *bc_ptr, uint8_t N_dof)
{
	void *bc = allocate_bc_atom(N_dof);
	bc_atom_copy(bc, bc_ptr, N_dof);
	return bc;	
}

void bc_atom_destroy(void *bc_ptr)
{
	bc_atom_finish(bc_ptr);
	nb_free_mem(bc_ptr);
}

void bc_atom_clear(void *bc_ptr, uint8_t N_dof)
{
	bc_atom_t *bc = bc_ptr;
	bc->id = 0;
	memset(bc->mask, 0, N_dof * sizeof(*(bc->mask)));
	memset(bc->val, 0, N_dof * sizeof(*(bc->val)));
}

void bc_atom_set_data(bc_atom_t *bc, uint32_t elem_id, uint8_t N_dof,
		      const bool dof_mask[], const double val[],
		      void (*fval)(const double *x, double t, double *out))
{
	bc->id = elem_id;
	memcpy(bc->mask, dof_mask, N_dof * sizeof(*(bc->mask)));
	if (NULL != fval)
		bc->fval = fval;
	else
		memcpy(bc->val, val, N_dof * sizeof(*(bc->val)));

}
