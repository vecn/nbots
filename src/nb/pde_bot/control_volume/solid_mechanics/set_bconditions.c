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

static uint32_t **malloc_elem_adj(const nb_partition_t *part);
static void free_elem_adj(uint32_t **elem_adj);
static void set_neumann_sgm(const nb_partition_t *part,
			    uint32_t **elem_adj, double* F, 
			    const nb_bcond_t *const bcond, 
			    double factor);
static void set_neumann_sgm_function(const nb_partition_t *part,
				     uint32_t **elem_adj,
				     double* F, uint8_t N_dof,
				     const nb_bcond_iter_t *const iter,
				     double factor);
static void set_neumann_sgm_integrated(const nb_partition_t *part,
				       uint32_t **elem_adj,
				       double* F, uint8_t N_dof,
				       const nb_bcond_iter_t *const iter,
				       double factor);
static void set_neumann(const nb_partition_t *part,
			uint8_t N_dof, double* F, double factor,
			double val[2], bool mask[2],
			uint32_t elem_id);
static void set_neumann_vtx(const nb_partition_t *part,
			    uint32_t **elem_adj, double* F, 
			    const nb_bcond_t *const bcond, 
			    double factor);
static void set_neumann_subsgm_adj_to_node(const nb_partition_t *part,
					   uint32_t **elem_adj,
					   double* F, uint8_t N_dof,
					   double factor, double val[2],
					   bool mask[2], uint32_t sgm_id,
					   uint32_t subsgm_id);
static void get_subsgm_adj_to_node(const nb_partition_t *part,
				   uint32_t node_id, uint32_t subsgm_id[4]);
static void set_dirichlet_sgm(const nb_partition_t *part,
			      uint32_t **elem_adj,
			      const nb_bcond_t *const bcond, 
			      double factor,
			      nb_bcond_t *numeric_bcond);
static void set_dirichlet_vtx(const nb_partition_t *part,
			      uint32_t **elem_adj,
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

void nb_cvfa_set_bconditions(const nb_partition_t *part,
			     vcn_sparse_t* K, double* F,
			     const nb_bcond_t *bcond,
			     double factor)
{
	uint16_t bcond_size = nb_bcond_get_memsize(2);
	nb_bcond_t *numeric_bcond = NB_SOFT_MALLOC(bcond_size);
	nb_bcond_init(numeric_bcond, 2);

	uint32_t **elem_adj = malloc_elem_adj(part);
	nb_partition_insgm_get_elem_adj(part, elem_adj);

	set_neumann_sgm(part, elem_adj, F, bcond, factor);
	set_neumann_vtx(part, elem_adj, F, bcond, factor);
	set_dirichlet_sgm(part, elem_adj, bcond, factor, numeric_bcond);
	set_dirichlet_vtx(part, elem_adj, bcond, factor, numeric_bcond);

	set_numeric_bconditions(K, F, part, numeric_bcond);

	free_elem_adj(elem_adj);
	NB_SOFT_FREE(bcond_size, numeric_bcond);
}

static uint32_t **malloc_elem_adj(const nb_partition_t *part)
{
	uint32_t N_insgm = nb_partition_get_N_insgm(part);
	uint32_t memsize = N_insgm * sizeof(uint32_t*);
	for (uint32_t i = 0; i < N_insgm; i++) {
		uint16_t N_subsgm = nb_partition_insgm_get_N_subsgm(part, i);
		memsize += N_subsgm * sizeof(uint32_t);
	}

	char *memblock = malloc(memsize);

	uint32_t **elem_adj = (void*) memblock;
	memblock += N_insgm * sizeof(uint32_t*);
	for (uint32_t i = 0; i < N_insgm; i++) {
		elem_adj[i] = (void*) memblock;
		uint16_t N_subsgm = nb_partition_insgm_get_N_subsgm(part, i);
		memblock += N_subsgm * sizeof(uint32_t);
	}

	return elem_adj;
}

static inline void free_elem_adj(uint32_t **elem_adj)
{
	free(elem_adj);
}

static void set_neumann_sgm(const nb_partition_t *part,
			    uint32_t **elem_adj,  double* F, 
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
			set_neumann_sgm_function(part, elem_adj, F,
						 N_dof, iter, factor);
		else
			set_neumann_sgm_integrated(part, elem_adj, F,
						   N_dof, iter, factor);
	}
	nb_bcond_iter_finish(iter);
}

static void set_neumann_sgm_function(const nb_partition_t *part,
				     uint32_t **elem_adj,
				     double* F, uint8_t N_dof,
				     const nb_bcond_iter_t *const iter,
				     double factor)
{
	uint32_t sgm_id = nb_bcond_iter_get_id(iter);

	uint32_t v1_id = nb_partition_get_node_x_insgm(part, sgm_id, 0);
	double x[2];
	x[0] = nb_partition_get_x_node(part, v1_id);
	x[1] = nb_partition_get_y_node(part, v1_id);
	double *val1 = alloca(N_dof * sizeof(double));
	nb_bcond_iter_get_val(iter, N_dof, x, 0, val1);

	double *val2 = alloca(N_dof * sizeof(double));

	uint32_t N = nb_partition_insgm_get_N_subsgm(part, sgm_id);
	for (uint32_t i = 0; i < N; i++) {
		double subsgm_length =
			nb_partition_insgm_subsgm_get_length(part, sgm_id, i);
		
		uint32_t v2_id = nb_partition_get_node_x_insgm(part, sgm_id,
							       i + 1);
		x[0] = nb_partition_get_x_node(part, v2_id);
		x[1] = nb_partition_get_y_node(part, v2_id);
		nb_bcond_iter_get_val(iter, N_dof, x, 0, val2);

		double val[2];
		val[0] = 0.5 * (val1[0] + val2[0]) * subsgm_length;
		val[1] = 0.5 * (val1[1] + val2[1]) * subsgm_length;

		bool mask[2] = {nb_bcond_iter_get_mask(iter, 0),
				nb_bcond_iter_get_mask(iter, 1)};

		uint32_t elem_id = elem_adj[sgm_id][i];
		set_neumann(part, N_dof, F, factor, val, mask, elem_id);

		v1_id = v2_id;
		memcpy(val1, val2, N_dof * sizeof(double));
	}
}

static void set_neumann_sgm_integrated(const nb_partition_t *part,
				       uint32_t **elem_adj,
				       double* F, uint8_t N_dof,
				       const nb_bcond_iter_t *const iter,
				       double factor)
{
	uint32_t model_id = nb_bcond_iter_get_id(iter);
	double sgm_length = nb_partition_insgm_get_length(part, model_id);

	uint32_t N = nb_partition_insgm_get_N_subsgm(part, model_id);
	for (uint32_t i = 0; i < N; i++) {
		double subsgm_length =
			nb_partition_insgm_subsgm_get_length(part, model_id, i);
		double w = subsgm_length / sgm_length;

		uint32_t elem_id = elem_adj[model_id][i];

		double x_dummy[2] = {0, 0};
		double *val = alloca(N_dof * sizeof(double));
		nb_bcond_iter_get_val(iter, N_dof, x_dummy, 0, val);

		bool mask[2] = {nb_bcond_iter_get_mask(iter, 0),
				nb_bcond_iter_get_mask(iter, 1)};
		set_neumann(part, N_dof, F, factor * w, val, mask, elem_id);
	}
}

static void set_neumann(const nb_partition_t *part,
			uint8_t N_dof, double* F, double factor,
			double val[2], bool mask[2],
			uint32_t elem_id)
{	
	for (uint8_t j = 0; j < N_dof; j++) {
		if (mask[j]) {
			uint32_t mtx_id = elem_id * N_dof + j;
			F[mtx_id] += factor * val[j];
		}
	}
}

static void set_neumann_vtx(const nb_partition_t *part,
			    uint32_t **elem_adj, double* F, 
			    const nb_bcond_t *const bcond, 
			    double factor)
{
	uint8_t N_dof = nb_bcond_get_N_dof(bcond);
	uint16_t size = nb_bcond_iter_get_memsize();
	nb_bcond_iter_t *iter = alloca(size);
	nb_bcond_iter_init(iter);
	nb_bcond_iter_set_conditions(iter, bcond, NB_NEUMANN,
				     NB_BC_ON_POINT);
	while (nb_bcond_iter_has_more(iter)) {
		nb_bcond_iter_go_next(iter);

		uint32_t node_id = nb_bcond_iter_get_id(iter);
		uint32_t subsgm_id[4];
		get_subsgm_adj_to_node(part, node_id, subsgm_id);

		double x[2];
		x[0] = nb_partition_get_x_node(part, node_id);
		x[1] = nb_partition_get_y_node(part, node_id);
		double *val = alloca(N_dof * sizeof(double));
		nb_bcond_iter_get_val(iter, N_dof, x, 0, val);

		bool mask[2] = {nb_bcond_iter_get_mask(iter, 0),
				nb_bcond_iter_get_mask(iter, 1)};

		set_neumann_subsgm_adj_to_node(part, elem_adj, F, N_dof,
					       factor, val, mask,
					       subsgm_id[0], subsgm_id[1]);

		set_neumann_subsgm_adj_to_node(part, elem_adj, F, N_dof,
					       factor, val, mask,
					       subsgm_id[2], subsgm_id[3]);
	}
	nb_bcond_iter_finish(iter);
}

static void set_neumann_subsgm_adj_to_node(const nb_partition_t *part,
					   uint32_t **elem_adj,
					   double* F, uint8_t N_dof,
					   double factor, double val[2],
					   bool mask[2], uint32_t sgm_id,
					   uint32_t subsgm_id)
{
	double subsgm_length =
		nb_partition_insgm_subsgm_get_length(part, sgm_id, subsgm_id);
	uint32_t elem_id = elem_adj[sgm_id][subsgm_id];
	double integral_factor = 0.5 * subsgm_length;
	set_neumann(part, N_dof, F, factor * integral_factor,
		    val, mask, elem_id);
}

static void get_subsgm_adj_to_node(const nb_partition_t *part,
				   uint32_t node_id, uint32_t subsgm_id[4])
{
	uint32_t N_sgm = nb_partition_get_N_insgm(part);
	for (uint32_t i = 0; i < N_sgm; i++) {
		uint32_t N = nb_partition_get_N_nodes_x_insgm(part, i);
		for (uint16_t j = 0; j < N; j++) {
			uint32_t sid = nb_partition_get_node_x_insgm(part, i, j);
			if (sid == node_id) {
				if (j > 0) {
					subsgm_id[i * 2] = i;
					subsgm_id[i*2+1] = j - 1;
					i ++;
				}

				if (j < N - 1) {
					subsgm_id[i * 2] = i;
					subsgm_id[i*2+1] = j;
					i ++;
				}
				
				if (i < 2)
					break;
				else
					goto EXIT;
			}
		}
	}
EXIT:
	return;
}

static void set_dirichlet_sgm(const nb_partition_t *part,
			      uint32_t **elem_adj,
			      const nb_bcond_t *const bcond, 
			      double factor,
			      nb_bcond_t *numeric_bcond)
{
	uint8_t N_dof = nb_bcond_get_N_dof(bcond);
	uint16_t size = nb_bcond_iter_get_memsize();
	nb_bcond_iter_t *iter = alloca(size);
	nb_bcond_iter_init(iter);
	nb_bcond_iter_set_conditions(iter, bcond, NB_DIRICHLET,
				     NB_BC_ON_SEGMENT);
	while (nb_bcond_iter_has_more(iter)) {
		nb_bcond_iter_go_next(iter);
		uint32_t sgm_id = nb_bcond_iter_get_id(iter);
		uint32_t N = nb_partition_insgm_get_N_subsgm(part, sgm_id);
		for (uint32_t i = 0; i < N; i++) {
			uint32_t elem_id = elem_adj[sgm_id][i];
			
			bool mask[2] = {nb_bcond_iter_get_mask(iter, 0),
					nb_bcond_iter_get_mask(iter, 1)};
			
			double x[2] = {nb_partition_get_x_elem(part, elem_id),
				       nb_partition_get_y_elem(part, elem_id)};
			double *val = alloca(N_dof * sizeof(double));
			nb_bcond_iter_get_val(iter, N_dof, x, 0, val);
			val[0] *= factor;
			val[1] *= factor;

			nb_bcond_push(numeric_bcond, NB_DIRICHLET,
				      NB_BC_ON_POINT, elem_id, mask, val);
		}
	}
	nb_bcond_iter_finish(iter);	
}

static void set_dirichlet_vtx(const nb_partition_t *part,
			      uint32_t **elem_adj,
			      const nb_bcond_t *const bcond, 
			      double factor,
			      nb_bcond_t *numeric_bcond)
{
	uint8_t N_dof = nb_bcond_get_N_dof(bcond);
	uint16_t size = nb_bcond_iter_get_memsize();
	nb_bcond_iter_t *iter = alloca(size);
	nb_bcond_iter_init(iter);
	nb_bcond_iter_set_conditions(iter, bcond, NB_DIRICHLET,
				     NB_BC_ON_POINT);
	while (nb_bcond_iter_has_more(iter)) {
		nb_bcond_iter_go_next(iter);
		uint32_t node_id = nb_bcond_iter_get_id(iter);
		uint32_t subsgm_id[4];
		get_subsgm_adj_to_node(part, node_id, subsgm_id);
			
		bool mask[2] = {nb_bcond_iter_get_mask(iter, 0),
				nb_bcond_iter_get_mask(iter, 1)};
			
		double x[2] = {nb_partition_get_x_node(part, node_id),
			       nb_partition_get_y_node(part, node_id)};
		double *val = alloca(N_dof * sizeof(double));
		nb_bcond_iter_get_val(iter, N_dof, x, 0, val);
		val[0] *= factor;
		val[1] *= factor;

		uint8_t N_subsgm_adj_to_node = 2;
		for (uint8_t i = 0; i < N_subsgm_adj_to_node; i++) {
			uint32_t elem_id = 
				elem_adj[subsgm_id[i*2]][subsgm_id[i*2+1]];
			nb_bcond_push(numeric_bcond, NB_DIRICHLET,
				      NB_BC_ON_POINT, elem_id, mask, val);
		}
	}
	nb_bcond_iter_finish(iter);
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