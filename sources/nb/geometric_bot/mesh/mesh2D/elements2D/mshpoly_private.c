#include <stdint.h>

#include "mshpoly_struct.h"
typedef struct nb_mshpoly_s nb_mshpoly_t;

void nb_mshpoly_node_move_x(void *msh, uint32_t i, double x)
{
	nb_mshpoly_t *mshpoly = msh;
	mshpoly->nod[i * 2] += x;
}

void nb_mshpoly_node_move_y(void *msh, uint32_t i, double y)
{
	nb_mshpoly_t *mshpoly = msh;
	mshpoly->nod[i*2+1] += y;
}

void nb_mshpoly_elem_move_x(void *msh, uint32_t i, double x)
{
	nb_mshpoly_t *mshpoly = msh;
	mshpoly->cen[i * 2] += x;
}

void nb_mshpoly_elem_move_y(void *msh, uint32_t i, double y)
{
	nb_mshpoly_t *mshpoly = msh;
	mshpoly->cen[i*2+1] += y;
}
