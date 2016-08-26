#include <stdlib.h>
#include <stdint.h>

#include "nb/geometric_bot/mesh/partition.h"
#include "partition_struct.h"

static void assemble_system(const nb_partition_t *part,
			    uint8_t N_comp, double *M, double *F,
			    const double *elem_values);
static void solve_system(const nb_partition_t *part,
			 uint8_t N_comp, double *M, double *F,
			 double *nodal_values);

void nb_partition_extrapolate_elems_to_nodes(const nb_partition_t *part,
					     uint8_t N_comp,
					     const double *elem_values,
					     double *nodal_values)
{
	uint32_t N_nodes = part->get_N_nodes(part->msh);

	char *memblock = malloc(2 * N_nodes * N_comp * sizeof(double));

	double *M = (void*) memblock;
	double *F = (void*) (memblock + N_nodes * N_comp * sizeof(double));
	memset(M, 0, N_nodes * N_comp * sizeof(*M));
	memset(F, 0, N_nodes * N_comp * sizeof(*M));
	
	assemble_system(part, N_comp, M, F, elem_values);

	solve_system(part, N_comp, M, F, nodal_values);

	free(M);
}

static void assemble_system(const nb_partition_t *part,
			    uint8_t N_comp, double *M, double *F,
			    const double *elem_values)
{
	uint32_t N_elems = part->get_N_elems(part->msh);
	/* AQUI VOY: distort y draw func to partition implementations */
}

static void solve_system(const nb_partition_t *part,
			 uint8_t N_comp, double *M, double *F,
			 double *nodal_values)
{
	uint32_t N_nodes = part->get_N_nodes(part->msh);
	uint32_t N = N_comp * N_nodes;
	for (uint32_t i = 0; i < N; i++)
		nodal_values[i] = F[i] / M[i];
}

void nb_partition_extrapolate_elem_points_to_nodes
				(const nb_partition_t *part,
				 uint8_t N_comp,
				 const double *elem_values,
				 double *nodal_values,
				 void (*assemble_elem)(const nb_partition_t *part,
						       uint8_t N_comp,
						       const double *elem_values,
						       uint32_t elem_id,
						       double *F, double *M))
{
	uint32_t N_nodes = part->get_N_nodes(part->msh);

	char *memblock = malloc(2 * N_nodes * N_comp * sizeof(double));

	double *M = (void*) memblock;
	double *F = (void*) (memblock + N_nodes * N_comp * sizeof(double));
	memset(M, 0, N_nodes * N_comp * sizeof(*M));
	memset(F, 0, N_nodes * N_comp * sizeof(*M));
	
	assemble_system(part, N_comp, M, F, elem_values, assemble_elem);

	solve_system(part, N_comp, M, F, nodal_values);

	free(M);
}
