#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <alloca.h>

#include "nb/eigen_bot.h"
#include "nb/geometric_bot.h"
#include "nb/pde_bot/boundary_conditions/bcond.h"
#include "nb/pde_bot/boundary_conditions/bcond_iter.h"

#include "set_bconditions.h"

static void set_numeric_bconditions(vcn_sparse_t *K, double *F,
				    const nb_partition_t *const part,
				    nb_bcond_t *bcond);

static void set_numeric_bcond_dirichlet(const nb_partition_t *part,
					vcn_sparse_t* K, uint8_t N_dof,
					double* F,
					nb_bcond_iter_t *iter,
					uint32_t elem_id);

void nb_cvfa_set_bcondtions(const nb_partition_t *part,
				   const nb_cond_t *bcond,
				   vcn_sparse_t* K, double* F,
				   double factor)
{
	uint16_t bcond_size = nb_bcond_get_memsize(2);
	nb_bcond_t *numeric_bcond = NB_SOFT_MALLOC(bcond_size);
	nb_bcond_init(numeric_bcond, 2);
	/* PENDING */

	set_numeric_bconditions(K, F, part, numeric_bcond);

	NB_SOFT_FREE(bcond_size, numeric_bcond);
}

static void set_numeric_bconditions(vcn_sparse_t *K, double *F,
				    const nb_partition_t *const part,
				    nb_bcond_t *bcond)
{
	uint8_t N_dof = nb_bcond_get_N_dof(bcond);
	uint16_t size = nb_bcond_iter_get_memsize();
	nb_bcond_iter_t *iter = alloca(size);
	nb_bcond_iter_init(iter);
	nb_bcond_iter_set_conditions(iter, bcond,
				     NB_DIRICHLET,
				     NB_BC_ON_POINT);
	while (nb_bcond_iter_has_more(iter)) {
		nb_bcond_iter_go_next(iter);
		uint32_t elem_id = nb_bcond_iter_get_id(iter);
		set_numeric_bcond_dirichlet(part, K, N_dof,
					    F, iter, elem_id);
	}
	nb_bcond_iter_finish(iter);
}

static void set_numeric_bcond_dirichlet(const nb_partition_t *part,
					vcn_sparse_t* K, uint8_t N_dof,
					double* F,
					nb_bcond_iter_t *iter,
					uint32_t elem_id)
{
	double x[2];
	x[0] = nb_partition_get_x_elem(part, elem_id);
	x[1] = nb_partition_get_y_elem(part, elem_id);

	double *val = alloca(N_dof * sizeof(*val));
	nb_bcond_iter_get_val(iter, N_dof, x, 0, val);
	for (uint8_t j = 0; j < N_dof; j++) {
		bool mask = nb_bcond_iter_get_mask(iter, j);
		if (mask) {
			uint32_t mtx_id = elem_id * N_dof + j;
			vcn_sparse_set_Dirichlet_condition(K, F, mtx_id,
							   val[j]);
		}
	}
}
