#include <stdint.h>

#include "nb/geometric_bot.h"

#include "msh3trg_struct.h"

typedef struct nb_msh3trg_s nb_msh3trg_t;

void nb_msh3trg_node_move_x(void *msh, uint32_t i, double x)
{
	nb_msh3trg_t *msh3trg = msh;
	msh3trg->nod[i * 2] += x;
}

void nb_msh3trg_node_move_y(void *msh, uint32_t i, double y)
{
	nb_msh3trg_t *msh3trg = msh;
	msh3trg->nod[i*2+1] += y;
}

void nb_msh3trg_elem_move_x(void *msh, uint32_t i, double x)
{
	;/* Do nothing */
}

void nb_msh3trg_elem_move_y(void *msh, uint32_t i, double y)
{
	;/* Do nothing */
}

uint32_t nb_msh3trg_get_N_total_adj(const void *msh)
{
	uint32_t N_elems = nb_msh3trg_get_N_elems(msh);
	return N_elems * 3;
}
