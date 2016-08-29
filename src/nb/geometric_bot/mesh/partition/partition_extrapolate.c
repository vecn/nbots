#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "nb/memory_bot.h"
#include "nb/geometric_bot/mesh/partition.h"
#include "partition_struct.h"

#define POW2(a) ((a)*(a))

static void extrapolate_to_nodes
		(const nb_partition_t *part,
		 uint8_t N_comp,
		 const double *elem_values,
		 double *nodal_values,
		 void (*assemble_elem)(const nb_partition_t *part,
				       const double *elem_values,
				       uint32_t elem_id,
				       uint8_t N_comp,
				       double *M, double *F));
static void assemble_system(const nb_partition_t *part,
			    const double *elem_values,
			    uint8_t N_comp, double *M, double *F,
			    void (*assemble_elem)(const nb_partition_t *part,
						  const double *elem_values,
						  uint32_t elem_id,
						  uint8_t N_comp,
						  double *M, double *F));
static void solve_system(const nb_partition_t *part,
			 uint8_t N_comp, double *M, double *F,
			 double *nodal_values);
static void assemble_elem_centroids(const nb_partition_t *part,
				    const double *elem_values,
				    uint32_t elem_id,
				    uint8_t N_comp,
				    double *M, double *F);
static void eval_shape_funcs(const nb_partition_t *part, uint32_t elem_id,
			     double *f);
static double monotonic_func(double x);

void nb_partition_extrapolate_elems_to_nodes(const nb_partition_t *part,
					     uint8_t N_comp,
					     const double *elem_values,
					     double *nodal_values)
{
	extrapolate_to_nodes(part, N_comp, elem_values, nodal_values,
			     assemble_elem_centroids);
}

static void extrapolate_to_nodes
		(const nb_partition_t *part,
		 uint8_t N_comp,
		 const double *elem_values,
		 double *nodal_values,
		 void (*assemble_elem)(const nb_partition_t *part,
				       const double *elem_values,
				       uint32_t elem_id,
				       uint8_t N_comp,
				       double *M, double *F))
{
	uint32_t N_nodes = part->get_N_nodes(part->msh);

	char *memblock = malloc((1 + N_comp) * N_nodes * sizeof(double));

	double *M = (void*) memblock;
	double *F = (void*) (memblock + N_nodes * sizeof(double));
	memset(M, 0, N_nodes * sizeof(*M));
	memset(F, 0, N_nodes * N_comp * sizeof(*F));
	
	assemble_system(part, elem_values, N_comp, M, F, assemble_elem);

	solve_system(part, N_comp, M, F, nodal_values);

	free(memblock);
}

static void assemble_system(const nb_partition_t *part,
			    const double *elem_values,
			    uint8_t N_comp, double *M, double *F,
			    void (*assemble_elem)(const nb_partition_t *part,
						  const double *elem_values,
						  uint32_t elem_id,
						  uint8_t N_comp,
						  double *M, double *F))
{
	uint32_t N_elems = part->get_N_elems(part->msh);
	for (uint32_t i = 0; i < N_elems; i++)
		assemble_elem(part, elem_values, i, N_comp, M, F);
}

static void solve_system(const nb_partition_t *part,
			 uint8_t N_comp, double *M, double *F,
			 double *nodal_values)
{
	uint32_t N_nodes = part->get_N_nodes(part->msh);
	for (uint32_t i = 0; i < N_nodes; i++) {
		for (uint8_t j = 0; j < N_comp; j++) {
			uint32_t id = i * N_comp + j;
			nodal_values[id] = F[id] / M[i];
		}
	}
}

static void assemble_elem_centroids(const nb_partition_t *part,
				    const double *elem_values,
				    uint32_t elem_id,
				    uint8_t N_comp,
				    double *M, double *F)
{
	uint16_t N_adj = part->elem_get_N_adj(part->msh, elem_id);
	uint16_t f_memsize = N_adj * sizeof(double);
	double *f = NB_SOFT_MALLOC(f_memsize);
	eval_shape_funcs(part, elem_id, f);

	double area = part->elem_get_area(part->msh, elem_id);

	for (uint32_t i = 0; i < N_adj; i++) {
		uint32_t node_id = part->elem_get_adj(part->msh, elem_id, i);
		/* Assemble mass matrix */
		for (uint32_t j = 0; j < N_adj; j++) {
			M[node_id] += area * f[i] * f[j];
		}
		/* Assemble independent vector for each component */
		for (uint8_t j = 0; j < N_comp; j++) {
			uint32_t nid = node_id * N_comp + j;
			uint32_t eid = elem_id * N_comp + j;
			F[nid] += area * elem_values[eid] * f[i];
		}
	}

	NB_SOFT_FREE(f_memsize, f);
}

static void eval_shape_funcs(const nb_partition_t *part, uint32_t elem_id,
			     double *f)
{
	uint16_t N_adj = part->elem_get_N_adj(part->msh, elem_id);
	uint16_t memsize = 2 * N_adj * sizeof(double);
	char *memblock = NB_SOFT_MALLOC(memsize);
	
	double *g = (void*) memblock;
	double *h = (void*) (memblock + N_adj * sizeof(double));
	
	double x = part->elem_get_x(part->msh, elem_id);
	double y = part->elem_get_y(part->msh, elem_id);
	for (uint16_t i = 0; i < N_adj; i++) {
		double node_id = part->elem_get_adj(part->msh, elem_id, i);
		double xn = part->node_get_x(part->msh, node_id);
		double yn = part->node_get_x(part->msh, node_id);
		double dist2 = POW2(x - xn) + POW2(y - yn);
		h[i] = monotonic_func(dist2);
	}

	for (uint16_t i = 0; i < N_adj; i++) {
		g[i] = 1;
		for (uint16_t j = 0; j < N_adj; j++) {
			if (i != j)
				g[i] *= h[j];
		}
	}

	double sum = 0;
	for (uint16_t i = 0; i < N_adj; i++)
		sum += g[i];

	for (uint16_t i = 0; i < N_adj; i++)
		f[i] = g[i] / sum;

	NB_SOFT_FREE(memsize, memblock);
}

static inline double monotonic_func(double x)
{
	return x;
}
