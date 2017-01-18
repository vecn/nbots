#include <stdint.h>

#include "nb/geometric_bot.h"

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

uint32_t nb_mshpoly_get_N_total_adj(const void *msh)
{
	uint32_t N_elems = nb_mshpoly_get_N_elems(msh);
	uint32_t N_total_adj = 0;
	for (uint32_t i = 0; i < N_elems; i++) {
		uint32_t N_adj = nb_mshpoly_elem_get_N_adj(msh, i);
		N_total_adj += N_adj;
	}
	return N_total_adj;
}
