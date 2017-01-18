#include <stdint.h>

#include "nb/geometric_bot.h"

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

uint32_t nb_mshquad_get_N_total_adj(const void *msh)
{
	uint32_t N_elems = nb_mshquad_get_N_elems(msh);
	uint32_t N_total_adj = 0;
	for (uint32_t i = 0; i < N_elems; i++) {
		if (nb_mshquad_elem_is_quad(msh, i))
			N_total_adj += 4;
		else
			N_total_adj += 3;
	}
	return N_total_adj;
}
