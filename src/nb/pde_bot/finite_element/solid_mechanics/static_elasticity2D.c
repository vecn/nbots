#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "nb/memory_bot.h"
#include "nb/math_bot.h"
#include "nb/cfreader_cat.h"
#include "nb/solver_bot.h"
#include "nb/container_bot.h"
#include "nb/graph_bot.h"
#include "nb/geometric_bot.h"
#include "nb/pde_bot/material.h"
#include "nb/pde_bot/common_solid_mechanics/analysis2D.h"
#include "nb/pde_bot/common_solid_mechanics/formulas.h"
#include "nb/pde_bot/boundary_conditions/bcond.h"
#include "nb/pde_bot/boundary_conditions/bcond_iter.h"
#include "nb/pde_bot/finite_element/element.h"
#include "nb/pde_bot/finite_element/gaussp_to_nodes.h"
#include "nb/pde_bot/finite_element/solid_mechanics/static_elasticity2D.h"

#include "set_bconditions.h"
#include "pipeline.h"

#define POW2(a) ((a)*(a))

static int solver(const nb_sparse_t *const A,
		  const double *const b, double* x);

int nb_fem_compute_2D_Solid_Mechanics
			(const nb_partition_t *const part,
			 const nb_fem_elem_t *const elemtype,
			 const nb_material_t *const material,
			 const nb_bcond_t *const bcond,
			 bool enable_self_weight,
			 double gravity[2],
			 nb_analysis2D_t analysis2D,
			 nb_analysis2D_params *params2D,
			 const bool *elements_enabled, /* NULL to enable all */
			 double *displacement, /* Output */
			 double *strain        /* Output */)
{
	int status = 0;
	nb_graph_t *graph = nb_allocate_mem(nb_graph_get_memsize());
	nb_graph_init(graph);
	nb_partition_load_graph(part, graph, NB_NODES_LINKED_BY_ELEMS);
	nb_sparse_t *K = nb_sparse_create(graph, NULL, 2);
	nb_graph_finish(graph);

	uint32_t N_nod = nb_partition_get_N_nodes(part);
	uint32_t F_memsize = 2 * N_nod * sizeof(double);
	double* F = nb_soft_allocate_mem(F_memsize);
	memset(F, 0, F_memsize);

	int status_assemble =
		pipeline_assemble_system(K, NULL, F, part, elemtype, material,
					 enable_self_weight, gravity,
					 analysis2D, params2D,
					 elements_enabled);
	if (0 != status_assemble) {
		status = 1;
		goto CLEANUP_LINEAR_SYSTEM;
	}

	nb_fem_set_bconditions(part, K, F, bcond, 1.0);

  
	int solver_status = solver(K, F, displacement);
	if (0 != solver_status) {
		status = 2;
		goto CLEANUP_LINEAR_SYSTEM;
	}

	pipeline_compute_strain(strain, part, displacement, elemtype);

CLEANUP_LINEAR_SYSTEM:
	nb_sparse_destroy(K);
	nb_soft_free_mem(F_memsize, F);
	return status;
}

static int solver(const nb_sparse_t *const A,
		  const double *const b, double* x)
{
	uint32_t N = nb_sparse_get_size(A);
	memset(x, 0, N * sizeof(*x));
	int status = nb_sparse_solve_CG_precond_Jacobi(A, b, x, N,
							1e-8, NULL,
							NULL, 1);
	int out;
	if (0 == status || 1 == status)
		out = 0;
	else
		out = 1; /* Tolerance not reached in CG Jacobi */
	return out;
}

void nb_fem_compute_stress_from_strain
			(uint32_t N_elements,
			 const nb_fem_elem_t *const elem,
			 const nb_material_t *const material,
			 nb_analysis2D_t analysis2D,
			 double* strain,
			 const bool* elements_enabled /* NULL to enable all */,
			 double* stress /* Output */)
{
	/* Compute stress from element strain */
	uint32_t omp_parallel_threads = 1;
#pragma omp parallel for num_threads(omp_parallel_threads) schedule(guided)
	for (uint32_t i = 0; i < N_elements; i++) {
		double D[4] = {1e-6, 1e-6, 1e-6, 1e-6};
		if (pipeline_elem_is_enabled(elements_enabled, i))
			nb_pde_get_constitutive_matrix(D, material,
						       analysis2D);

		uint8_t N_gp = nb_fem_elem_get_N_gpoints(elem);
		for (int j = 0; j < N_gp; j++) {
			uint32_t id = i * N_gp + j;
			stress[id * 3] = strain[id * 3] * D[0] +
				strain[id*3+1] * D[1];
			stress[id*3+1] = strain[id * 3] * D[1] +
				strain[id*3+1] * D[2];
			stress[id*3+2] = strain[id*3+2] * D[3];
		}
	}
}
