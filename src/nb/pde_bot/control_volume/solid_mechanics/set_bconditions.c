#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <alloca.h>

#include "nb/memory_bot.h"
#include "nb/eigen_bot.h"
#include "nb/geometric_bot.h"
#include "nb/pde_bot/boundary_conditions/bcond.h"
#include "nb/pde_bot/boundary_conditions/bcond_iter.h"

#include "set_bconditions.h"

static void set_neumann_sgm(const nb_partition_t *part,  double* F, 
			    const nb_bcond_t *const bcond, 
			    double factor);
static void set_neumann_sgm_function(const nb_partition_t *part,
				     double* F, uint8_t N_dof,
				     const nb_bcond_iter_t *const iter,
				     double factor);
static void set_neumann_sgm_integrated(const nb_partition_t *part,
				       double* F, uint8_t N_dof,
				       const nb_bcond_iter_t *const iter,
				       double factor);
static double nb_partition_insgm_get_length(const nb_partition_t *part,
				   uint32_t sgm_id);
static double nb_partition_insgm_subsgm_get_length(const nb_partition_t *part,
				      uint32_t sgm_id, uint32_t subsgm_id);
static void set_neumann_vtx(const nb_partition_t *part,
			    double* F, 
			    const nb_bcond_t *const bcond, 
			    double factor);
static void set_dirichlet_sgm(const nb_partition_t *part,
			      double* F, 
			      const nb_bcond_t *const bcond, 
			      double factor,
			      nb_bcond_t *numeric_bcond);
static void set_dirichlet_vtx(const nb_partition_t *part,
			      double* F, 
			      const nb_bcond_t *const bcond, 
			      double factor,
			      nb_bcond_t *numeric_bcond);
static void set_numeric_bconditions(vcn_sparse_t *K, double *F,
				    const nb_partition_t *const part,
				    nb_bcond_t *bcond);

static void set_numeric_bcond_dirichlet(const nb_partition_t *part,
					vcn_sparse_t* K, uint8_t N_dof,
					double* F,
					nb_bcond_iter_t *iter,
					uint32_t elem_id);

void nb_cvfa_set_bcondtions(const nb_partition_t *part,
			    const nb_bcond_t *bcond,
			    vcn_sparse_t* K, double* F,
			    double factor)
{
	uint16_t bcond_size = nb_bcond_get_memsize(2);
	nb_bcond_t *numeric_bcond = NB_SOFT_MALLOC(bcond_size);
	nb_bcond_init(numeric_bcond, 2);

	set_neumann_sgm(part, F, bcond, factor);
	set_neumann_vtx(part, F, bcond, factor);
	set_dirichlet_sgm(part, F, bcond, factor, numeric_bcond);
	set_dirichlet_vtx(part, F, bcond, factor, numeric_bcond);

	set_numeric_bconditions(K, F, part, numeric_bcond);

	NB_SOFT_FREE(bcond_size, numeric_bcond);
}

static void set_neumann_sgm(const nb_partition_t *part,  double* F, 
			    const nb_bcond_t *const bcond, 
			    double factor)
{
	uint8_t N_dof = nb_bcond_get_N_dof(bcond);
	uint16_t size = nb_bcond_iter_get_memsize();
	nb_bcond_iter_t *iter = alloca(size);
	nb_bcond_iter_init(iter);
	nb_bcond_iter_set_conditions(iter, bcond, NB_NEUMANN,
				     NB_BC_ON_SEGMENT);
	while (nb_bcond_iter_has_more(iter)) {
		nb_bcond_iter_go_next(iter);
		
		if (nb_bcond_iter_val_is_function(iter))
			set_neumann_sgm_function(part, F, N_dof,
						 iter, factor);
		else
			set_neumann_sgm_integrated(part, F, N_dof,
						   iter, factor);
	}
	nb_bcond_iter_finish(iter);
}

static void set_neumann_sgm_function(const nb_partition_t *part,
				     double* F, uint8_t N_dof,
				     const nb_bcond_iter_t *const iter,
				     double factor)
{
	/* PENDING */
	uint32_t model_id = nb_bcond_iter_get_id(iter);

	uint32_t v1_id = nb_partition_get_node_x_insgm(part, model_id, 0);
	double x[2];
	x[0] = nb_partition_get_x_node(part, v1_id);
	x[1] = nb_partition_get_y_node(part, v1_id);
	double *val1 = alloca(N_dof * sizeof(double));
	nb_bcond_iter_get_val(iter, N_dof, x, 0, val1);

	double *val2 = alloca(N_dof * sizeof(double));

	uint32_t N = nb_partition_get_N_nodes_x_insgm(part, model_id);
	for (uint32_t i = 0; i < N - 1; i++) {
		double subsgm_length =
			nb_partition_insgm_subsgm_get_length(part, model_id, i);
		
		uint32_t v2_id = nb_partition_get_node_x_insgm(part, model_id,
							       i + 1);
		x[0] = nb_partition_get_x_node(part, v2_id);
		x[1] = nb_partition_get_y_node(part, v2_id);
		nb_bcond_iter_get_val(iter, N_dof, x, 0, val2);

		for (uint8_t j = 0; j < N_dof; j++) {
			bool mask = nb_bcond_iter_get_mask(iter, j);
			if (mask) {
				double val = 0.5 * (val1[j] + val2[j]) *
					subsgm_length;

				uint32_t mtx_id1 = v1_id * N_dof + j;
				F[mtx_id1] += factor * val * 0.5;

				uint32_t mtx_id2 = v2_id * N_dof + j;
				F[mtx_id2] += factor * val * 0.5;
			}
		}

		v1_id = v2_id;
		memcpy(val1, val2, N_dof * sizeof(double));
	}
}

static void set_neumann_sgm_integrated(const nb_partition_t *part,
				       double* F, uint8_t N_dof,
				       const nb_bcond_iter_t *const iter,
				       double factor)
{
	/* PENDING */
	uint32_t model_id = nb_bcond_iter_get_id(iter);
	double sgm_length = nb_partition_insgm_get_length(part, model_id);

	uint32_t N = nb_partition_get_N_nodes_x_insgm(part, model_id);
	for (uint32_t i = 0; i < N - 1; i++) {
		double subsgm_length =
			nb_partition_insgm_subsgm_get_length(part, model_id, i);
		double w = subsgm_length / sgm_length;

		uint32_t v1_id = nb_partition_get_node_x_insgm(part,
							       model_id,
							       i);
		uint32_t v2_id = nb_partition_get_node_x_insgm(part,
							       model_id,
							       i + 1);

		set_neumann(part, N_dof, F, factor * w * 0.5, iter, v1_id);
		set_neumann(part, N_dof, F, factor * w * 0.5, iter, v2_id);
	}
}

static void set_neumann_vtx(const nb_partition_t *part,
			    double* F, 
			    const nb_bcond_t *const bcond, 
			    double factor)
{
	/* PENDING */

}

static void set_dirichlet_sgm(const nb_partition_t *part,
			      double* F, 
			      const nb_bcond_t *const bcond, 
			      double factor,
			      nb_bcond_t *numeric_bcond)
{
	/* PENDING */

}

static void set_dirichlet_vtx(const nb_partition_t *part,
			      double* F, 
			      const nb_bcond_t *const bcond, 
			      double factor,
			      nb_bcond_t *numeric_bcond)
{
	/* PENDING */

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
