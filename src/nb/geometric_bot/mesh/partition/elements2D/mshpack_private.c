#include <stdint.h>

#include "mshpack_struct.h"
typedef struct nb_mshpack_s nb_mshpack_t;

void nb_mshpack_node_move_x(void *msh, uint32_t i, double x)
{
	; /* Do nothing */
}

void nb_mshpack_node_move_y(void *msh, uint32_t i, double y)
{
	; /* Do nothing */
}

void nb_mshpack_elem_move_x(void *msh, uint32_t i, double x)
{
	nb_mshpack_t *mshpack = msh;
	mshpack->cen[i * 2] += x;
}

void nb_mshpack_elem_move_y(void *msh, uint32_t i, double y)
{
	nb_mshpack_t *mshpack = msh;
	mshpack->cen[i*2+1] += y;
}
