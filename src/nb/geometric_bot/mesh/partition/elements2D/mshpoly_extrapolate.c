#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "nb/memory_bot.h"
#include "nb/geometric_bot/mesh/partition/elements2D/mshpoly.h"

#include "mshpoly_struct.h"

#define POW2(a) ((a)*(a))

static void assemble_system(const void *part,
			    const double *elem_values,
			    uint8_t N_comp, double *M, double *F);
static void assemble_elem(const void *msh, const double *elem_values,
			  uint32_t elem_id, uint8_t N_comp,
			  double *M, double *F);
static void eval_shape_funcs(const void *msh, uint32_t elem_id, double *f);
static double monotonic_func(double x);
static void solve_system(const void *msh, uint8_t N_comp, 
			 double *M, double *F, double *nodal_values);
static double distort_using_nodal_field(nb_mshpoly_t *msh, double *disp,
					double max_disp);
static void elem_get_interpolation(const nb_mshpoly_t *msh, uint32_t elem_id,
				   const double *value, uint8_t N_comp,
				   double *out);
static double get_max_displacement(uint32_t N, double *disp);
static double distort_using_elem_field(nb_mshpoly_t *msh, double *disp,
				       double max_disp);

void nb_mshpoly_extrapolate_elems_to_nodes(const void *msh, uint8_t N_comp,
					   const double *elem_values,
					   double *nodal_values)
{
	uint32_t N_nodes = nb_mshpoly_get_N_nodes(msh);

	char *memblock = malloc((1 + N_comp) * N_nodes * sizeof(double));

	double *M = (void*) memblock;
	double *F = (void*) (memblock + N_nodes * sizeof(double));
	memset(M, 0, N_nodes * sizeof(*M));
	memset(F, 0, N_nodes * N_comp * sizeof(*F));
	
	assemble_system(msh, elem_values, N_comp, M, F);

	solve_system(msh, N_comp, M, F, nodal_values);

	free(memblock);
}

static void assemble_system(const void *msh,
			    const double *elem_values,
			    uint8_t N_comp, double *M, double *F)
{
	uint32_t N_elems = nb_mshpoly_get_N_elems(msh);
	for (uint32_t i = 0; i < N_elems; i++)
		assemble_elem(msh, elem_values, i, N_comp, M, F);
}

static void assemble_elem(const void *msh, const double *elem_values,
			  uint32_t elem_id, uint8_t N_comp,
			  double *M, double *F)
{
	uint16_t N_adj = nb_mshpoly_elem_get_N_adj(msh, elem_id);
	uint16_t f_memsize = N_adj * sizeof(double);
	double *f = NB_SOFT_MALLOC(f_memsize);
	eval_shape_funcs(msh, elem_id, f);

	double area = nb_mshpoly_elem_get_area(msh, elem_id);

	for (uint32_t i = 0; i < N_adj; i++) {
		uint32_t node_id = nb_mshpoly_elem_get_adj(msh, elem_id, i);
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

static void eval_shape_funcs(const void *msh, uint32_t elem_id,
			     double *f)
{
	uint16_t N_adj = nb_mshpoly_elem_get_N_adj(msh, elem_id);
	uint16_t memsize = 2 * N_adj * sizeof(double);
	char *memblock = NB_SOFT_MALLOC(memsize);
	
	double *g = (void*) memblock;
	double *h = (void*) (memblock + N_adj * sizeof(double));
	
	double x = nb_mshpoly_elem_get_x(msh, elem_id);
	double y = nb_mshpoly_elem_get_y(msh, elem_id);
	for (uint16_t i = 0; i < N_adj; i++) {
		double node_id = nb_mshpoly_elem_get_adj(msh, elem_id, i);
		double xn = nb_mshpoly_node_get_x(msh, node_id);
		double yn = nb_mshpoly_node_get_x(msh, node_id);
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

static void solve_system(const void *msh, uint8_t N_comp, 
			 double *M, double *F, double *nodal_values)
{
	uint32_t N_nodes = nb_msh3trg_get_N_nodes(msh);
	for (uint32_t i = 0; i < N_nodes; i++) {
		for (uint8_t j = 0; j < N_comp; j++) {
			uint32_t id = i * N_comp + j;
			nodal_values[id] = F[id] / M[i];
		}
	}
}


double nb_mshpoly_distort_with_field(void *msh,
				     nb_partition_entity field_entity,
				     double *disp,
				     double max_disp)
{
	double scale = 1.0;
	if (NB_NODE == field_entity)
		scale = distort_using_nodal_field(msh, disp, max_disp);
	else if (NB_ELEMENT == field_entity)
		scale = distort_using_elem_field(msh, disp, max_disp);
	return scale;
}

static double distort_using_nodal_field(nb_mshpoly_t *msh, double *disp,
					double max_disp)
{
	uint32_t N_nodes = nb_mshpoly_get_N_nodes(msh);
	double max_field_disp = get_max_displacement(N_nodes, disp);
	double scale = max_disp / max_field_disp;
	
	for (uint32_t i = 0; i < 2 * N_nodes; i++)
		msh->nod[i] += disp[i] * scale;


	uint32_t N_elems = nb_mshpoly_get_N_elems(msh);
	for (uint32_t i = 0; i < N_elems; i++) {
		double d[2];
		elem_get_interpolation(msh, i, disp, 2, d);
		msh->cen[i * 2] += d[0] * scale;
		msh->cen[i*2+1] += d[1] * scale;
	}

	return scale;
}

static void elem_get_interpolation(const nb_mshpoly_t *msh, uint32_t elem_id,
				   const double *value, uint8_t N_comp,
				   double *out)
{
	memset(out, 0, N_comp * sizeof(*out));
	
	uint32_t N_adj = nb_mshpoly_elem_get_N_adj(msh, elem_id);
	double *f = alloca(N_adj * sizeof(*f));
	eval_shape_funcs(msh, elem_id, f);

	for (uint32_t j = 0; j < N_adj; j++) {
		uint32_t id = nb_mshpoly_elem_get_adj(msh, elem_id, j);
		for (uint8_t c = 0; c < N_comp; c++)
			out[c] += value[id * N_comp + c] * f[j];
	}
}

static double get_max_displacement(uint32_t N, double *disp)
{
	double max = 0;
	for (uint32_t i = 0; i < N; i++) {
		double disp2 = POW2(disp[i * 2]) + POW2(disp[i*2+1]);
		if (disp2 > max)
			max = disp2;
	}
	return sqrt(max);
}

static double distort_using_elem_field(nb_mshpoly_t *msh, double *disp,
				       double max_disp)
{
	uint32_t N_nodes = nb_mshpoly_get_N_nodes(msh);
	uint32_t N_elems = nb_mshpoly_get_N_elems(msh);

	double max_field_disp = get_max_displacement(N_elems, disp);
	double scale = max_disp / max_field_disp;

	for (uint32_t i = 0; i < 2 * N_elems; i++)
		msh->cen[i] += disp[i] * scale;

	uint32_t memsize = 2 * N_nodes * sizeof(double);
	double *nodal_disp = NB_SOFT_MALLOC(memsize);
	nb_mshpoly_extrapolate_elems_to_nodes(msh, 2, disp, nodal_disp);
	
	for (uint32_t i = 0; i < 2 * N_nodes; i++)
		msh->nod[i] += nodal_disp[i] * scale;

	NB_SOFT_FREE(memsize, nodal_disp);
	return scale;
}
