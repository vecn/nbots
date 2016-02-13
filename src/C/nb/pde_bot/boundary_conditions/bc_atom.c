#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "bc_atom.h"

inline uint16_t bc_atom_get_memsize(uint8_t N_dof)
{
	uint16_t size = sizeof(*bc);
	uint16_t size_mask = sizeof(*(bc->mask));
	uint16_t size_val = sizeof(*(bc->val));
	return size + N_dof * (size_mask + size_val);
}

void* bc_atom_create(uint8_t N_dof)
{
	uint16_t size = sizeof(bc_atom_t);
	uint16_t size_mask = sizeof(*(bc->mask));
	uint16_t size_val = sizeof(*(bc->val));
	char *memblock = calloc(size + N_dof * (size_mask + size_val), 1);
	bc_atom_t *bc = memblock;
	bc->mask = memblock + size;
	bc->val = memblock + size + N_dof * size_mask;
}

inline void* bc_atom_clone(const void *const bc_ptr, uint8_t N_dof)
{
	const bc_atom_t *const bc = bc_ptr;
	uint16_t total_size = bc_atom_get_memsize(N_dof);
	char *bc_clone = calloc(total_size, 1);
	memcpy(bc_clone, (void*)bc, total_size);

}

inline void bc_atom_destroy(void *bc_ptr)
{
	free(bc);
}

inline void bc_atom_set_data(bc_atom_t *bc,uint32_t elem_id,
			     const bool dof_mask[],
			     const double value[])
{
	bc->id = elem_id;
	memcpy(bc->mask, dof_mask, bcond->N_dof * sizeof(*(bc->mask)));
	memcpy(bc->val, value, bcond->N_dof * sizeof(*(bc->val)));
}
