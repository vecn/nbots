#include <stdlib.h>
#include <stdint.h>

#include "nb/geometric_bot/mesh/partition.h"

void nb_partition_extrapolate_elems_to_nodes(nb_partition_t *part,
					     uint8_t N_comp,
					     const double *elem_values,
					     double *nodal_values)
{
	uint32_t N_elems = part->get_N_elems(part->msh);
	uint32_t N_nodes = part->get_N_nodes(part->msh);
	double *M = malloc(N_nod * N_comp * sizeof(*M));
	memset(M, 0, N_nod * N_comp * sizeof(*M));
/* AQUI VOY -> falta interfaz de partitions-> distort and draws functions
	free(M);
}
