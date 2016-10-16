#include <stdint.h>

#include "mshquad_struct.h"
typedef struct nb_mshquad_s nb_mshquad_t;

void nb_mshquad_node_move_x(void *msh, uint32_t i, double x)
{
	nb_mshquad_t *mshquad = msh;
	mshquad->nod[i * 2] += x;
}

void nb_mshquad_node_move_y(void *msh, uint32_t i, double y)
{
	nb_mshquad_t *mshquad = msh;
	mshquad->nod[i*2+1] += y;
}

void nb_mshquad_elem_move_x(void *msh, uint32_t i, double x)
{
	;/* Do nothing */
}

void nb_mshquad_elem_move_y(void *msh, uint32_t i, double y)
{
	;/* Do nothing */
}
