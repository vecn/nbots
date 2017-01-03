#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "nb/memory_bot.h"
#include "nb/solver_bot.h"
#include "nb/geometric_bot.h"
#include "nb/pde_bot/material.h"
#include "nb/pde_bot/common_solid_mechanics/analysis2D.h"
#include "nb/pde_bot/common_solid_mechanics/formulas.h"
#include "nb/pde_bot/boundary_conditions/bcond.h"
#include "nb/pde_bot/boundary_conditions/bcond_iter.h"

#include "set_bconditions.h"

#define POW2(a) ((a)*(a))

static uint32_t **allocate_elem_adj(const nb_mesh2D_t *mesh);
static void get_elem_adj(const nb_mesh2D_t *mesh, uint32_t **elem_adj);
static void check_elem_adj(const nb_mesh2D_t *mesh,
			   uint32_t **elem_adj, uint32_t elem_id);
static void check_boundary_face_adj(const nb_mesh2D_t *mesh,
				    uint32_t **elem_adj,
				    uint32_t elem_id, uint16_t face_id);
static bool face_is_the_same(uint32_t n1, uint32_t n2,
			     uint32_t s1, uint32_t s2);
static void free_elem_adj(uint32_t **elem_adj);
static void set_neumann_sgm(const nb_mesh2D_t *mesh,
			    const nb_material_t *material,
			    nb_analysis2D_t analysis2D,
			    uint32_t **elem_adj,
			    const nb_bcond_t *const bcond, 
			    double factor,
			    nb_bcond_t *numeric_bcond);
static void set_neumann_sgm_function(const nb_mesh2D_t *mesh,
				     const nb_material_t *material,
				     nb_analysis2D_t analysis2D,
				     uint32_t **elem_adj,
				     uint8_t N_dof,
				     const nb_bcond_iter_t *const iter,
				     double factor,
				     nb_bcond_t *numeric_bcond);
static void set_neumann_sgm_integrated(const nb_mesh2D_t *mesh,
				       uint32_t **elem_adj,
				       uint8_t N_dof,
				       const nb_bcond_iter_t *const iter,
				       double factor,
				       nb_bcond_t *numeric_bcond);
static void set_neumann_vtx(const nb_mesh2D_t *mesh,
			    uint32_t **elem_adj, 
			    const nb_bcond_t *const bcond, 
			    double factor,
			    nb_bcond_t *numeric_bcond);
static void set_neumann_subsgm_adj_to_node(const nb_mesh2D_t *mesh,
					   uint32_t **elem_adj,
					   uint8_t N_dof,
					   double factor, double val[2],
					   bool mask[2], uint32_t sgm_id,
					   uint32_t subsgm_id,
					   nb_bcond_t *numeric_bcond);
static void get_subsgm_adj_to_node(const nb_mesh2D_t *mesh,
				   uint32_t node_id, uint32_t subsgm_id[4]);
static void set_dirichlet_sgm(const nb_mesh2D_t *mesh,
			      uint32_t **elem_adj,
			      const nb_bcond_t *const bcond, 
			      double factor,
			      nb_bcond_t *numeric_bcond);
static void set_dirichlet_vtx(const nb_mesh2D_t *mesh,
			      uint32_t **elem_adj,
			      const nb_bcond_t *const bcond, 
			      double factor,
			      nb_bcond_t *numeric_bcond);
static void set_numeric_bcond_neumann(const nb_mesh2D_t *mesh,
				      uint32_t elem_id, 
				      uint8_t N_dof, double* F,
				      nb_bcond_iter_t *iter);
static void set_numeric_bcond_dirichlet(const nb_mesh2D_t *mesh,
					nb_sparse_t* K, uint8_t N_dof,
					double* F,
					nb_bcond_iter_t *iter,
					uint32_t elem_id);


void nb_cvfa_set_bconditions(const nb_mesh2D_t *mesh,
			     const nb_material_t *material,
			     nb_analysis2D_t analysis2D,
			     nb_sparse_t* K, double* F,
			     const nb_bcond_t *bcond,
			     double factor)
{
	uint16_t bcond_size = nb_bcond_get_memsize(2);
	nb_bcond_t *numeric_bcond = nb_soft_allocate_mem(bcond_size);
	nb_bcond_init(numeric_bcond, 2);

	nb_cvfa_get_numeric_bconditions(mesh, material, analysis2D, bcond,
					factor, numeric_bcond);
	nb_cvfa_set_numeric_bconditions(K, F, mesh, numeric_bcond);

	nb_bcond_finish(numeric_bcond);
	nb_soft_free_mem(bcond_size, numeric_bcond);
}

void nb_cvfa_get_numeric_bconditions(const nb_mesh2D_t *mesh,
				     const nb_material_t *material,
				     nb_analysis2D_t analysis2D,
				     const nb_bcond_t *const bcond,
				     double factor,
				     nb_bcond_t *numeric_bcond)
{
	uint32_t **elem_adj = allocate_elem_adj(mesh);
	get_elem_adj(mesh, elem_adj);

	set_neumann_sgm(mesh, material, analysis2D,
			elem_adj, bcond, factor, numeric_bcond);
	set_neumann_vtx(mesh, elem_adj, bcond, factor, numeric_bcond);
	set_dirichlet_sgm(mesh, elem_adj, bcond, factor, numeric_bcond);
	set_dirichlet_vtx(mesh, elem_adj, bcond, factor, numeric_bcond);

	free_elem_adj(elem_adj);
}

static uint32_t **allocate_elem_adj(const nb_mesh2D_t *mesh)
{
	uint32_t N_insgm = nb_mesh2D_get_N_insgm(mesh);
	uint32_t memsize = N_insgm * sizeof(uint32_t*);
	for (uint32_t i = 0; i < N_insgm; i++) {
		uint16_t N_subsgm = nb_mesh2D_insgm_get_N_subsgm(mesh, i);
		memsize += N_subsgm * sizeof(uint32_t);
	}

	char *memblock = nb_allocate_mem(memsize);

	uint32_t **elem_adj = (void*) memblock;
	memblock += N_insgm * sizeof(uint32_t*);
	for (uint32_t i = 0; i < N_insgm; i++) {
		elem_adj[i] = (void*) memblock;
		uint16_t N_subsgm = nb_mesh2D_insgm_get_N_subsgm(mesh, i);
		memblock += N_subsgm * sizeof(uint32_t);
	}

	return elem_adj;
}

static void get_elem_adj(const nb_mesh2D_t *mesh, uint32_t **elem_adj)
{
	uint32_t N_elems = nb_mesh2D_get_N_elems(mesh);
	for (uint32_t i = 0; i < N_elems; i++)
		check_elem_adj(mesh, elem_adj, i);
}

static void check_elem_adj(const nb_mesh2D_t *mesh,
			   uint32_t **elem_adj, uint32_t elem_id)
{
	uint16_t N_adj = nb_mesh2D_elem_get_N_adj(mesh, elem_id);
	for (uint16_t i = 0; i < N_adj; i++) {
		bool no_ngb = !nb_mesh2D_elem_has_ngb(mesh, elem_id, i);
		if (no_ngb)
			check_boundary_face_adj(mesh, elem_adj, elem_id, i);
	}
}

static void check_boundary_face_adj(const nb_mesh2D_t *mesh,
				    uint32_t **elem_adj,
				    uint32_t elem_id, uint16_t face_id)
{
	uint16_t N_adj = nb_mesh2D_elem_get_N_adj(mesh, elem_id);
	uint32_t n1 = nb_mesh2D_elem_get_adj(mesh, elem_id, face_id);
	uint32_t n2 = nb_mesh2D_elem_get_adj(mesh, elem_id,
						(face_id + 1) % N_adj);

	uint32_t N_insgm = nb_mesh2D_get_N_insgm(mesh);
	for (uint32_t i = 0; i < N_insgm; i++) {
		uint32_t N = nb_mesh2D_insgm_get_N_subsgm(mesh, i);
		for (uint32_t j = 0; j < N; j++) {
			uint32_t s1 = nb_mesh2D_insgm_get_node(mesh, i, j);
			uint32_t s2 = nb_mesh2D_insgm_get_node(mesh, i,
								  j + 1);
			if (face_is_the_same(n1, n2, s1, s2))
				elem_adj[i][j] = elem_id;
		}
	}
}

static inline bool face_is_the_same(uint32_t n1, uint32_t n2,
				    uint32_t s1, uint32_t s2)
{
	return ((n1 == s1) && (n2 == s2)) || ((n1 == s2) && (n2 == s1));
}

static inline void free_elem_adj(uint32_t **elem_adj)
{
	nb_free_mem(elem_adj);
}

static void set_neumann_sgm(const nb_mesh2D_t *mesh,
			    const nb_material_t *material,
			    nb_analysis2D_t analysis2D,
			    uint32_t **elem_adj,
			    const nb_bcond_t *const bcond, 
			    double factor,
			    nb_bcond_t *numeric_bcond)
{
	uint8_t N_dof = nb_bcond_get_N_dof(bcond);
	uint16_t size = nb_bcond_iter_get_memsize();
	nb_bcond_iter_t *iter = nb_allocate_on_stack(size);
	nb_bcond_iter_init(iter);
	nb_bcond_iter_set_conditions(iter, bcond, NB_NEUMANN,
				     NB_BC_ON_SEGMENT);
	while (nb_bcond_iter_has_more(iter)) {
		nb_bcond_iter_go_next(iter);
		
		if (nb_bcond_iter_val_is_function(iter))
			set_neumann_sgm_function(mesh, material, analysis2D,
						 elem_adj, N_dof, iter,
						 factor, numeric_bcond);
		else
			set_neumann_sgm_integrated(mesh, elem_adj, N_dof,
						   iter, factor, numeric_bcond);
	}
	nb_bcond_iter_finish(iter);
}

static void set_neumann_sgm_function(const nb_mesh2D_t *mesh,
				     const nb_material_t *material,
				     nb_analysis2D_t analysis2D,
				     uint32_t **elem_adj,
				     uint8_t N_dof,
				     const nb_bcond_iter_t *const iter,
				     double factor,
				     nb_bcond_t *numeric_bcond)
{
	uint32_t sgm_id = nb_bcond_iter_get_id(iter);

	uint32_t v1_id = nb_mesh2D_insgm_get_node(mesh, sgm_id, 0);
	double x[2];
	x[0] = nb_mesh2D_node_get_x(mesh, v1_id);
	x[1] = nb_mesh2D_node_get_y(mesh, v1_id);
	double *val1 = nb_allocate_on_stack(N_dof * sizeof(double));
	nb_bcond_iter_get_val(iter, N_dof, x, 0, val1);

	double *val2 = nb_allocate_on_stack(N_dof * sizeof(double));

	uint32_t N = nb_mesh2D_insgm_get_N_subsgm(mesh, sgm_id);
	for (uint32_t i = 0; i < N; i++) {
		double subsgm_length =
			nb_mesh2D_insgm_subsgm_get_length(mesh, sgm_id, i);
		
		uint32_t v2_id = nb_mesh2D_insgm_get_node(mesh, sgm_id,
							     i + 1);
		x[0] = nb_mesh2D_node_get_x(mesh, v2_id);
		x[1] = nb_mesh2D_node_get_y(mesh, v2_id);
		nb_bcond_iter_get_val(iter, N_dof, x, 0, val2);

		double val[2];
		val[0] = 0.5 * (val1[0] + val2[0]) * subsgm_length;
		val[1] = 0.5 * (val1[1] + val2[1]) * subsgm_length;

		double D[4];
		nb_pde_get_constitutive_matrix(D, material, analysis2D);

		bool mask[2] = {nb_bcond_iter_get_mask(iter, 0),
				nb_bcond_iter_get_mask(iter, 1)};

		uint32_t elem_id = elem_adj[sgm_id][i];

		val[0] *= factor;
		val[1] *= factor;
		nb_bcond_push(numeric_bcond, NB_NEUMANN, NB_BC_ON_POINT,
			      elem_id, mask, val);

		v1_id = v2_id;
		memcpy(val1, val2, N_dof * sizeof(*val1));
	}
}

static void set_neumann_sgm_integrated(const nb_mesh2D_t *mesh,
				       uint32_t **elem_adj,
				       uint8_t N_dof,
				       const nb_bcond_iter_t *const iter,
				       double factor,
				       nb_bcond_t *numeric_bcond)
{
	double *val = nb_allocate_on_stack(N_dof * sizeof(*val));

	uint32_t model_id = nb_bcond_iter_get_id(iter);
	uint32_t N = nb_mesh2D_insgm_get_N_subsgm(mesh, model_id);
	for (uint32_t i = 0; i < N; i++) {
		uint32_t elem_id = elem_adj[model_id][i];

		double x_dummy[2] = {0, 0};
		nb_bcond_iter_get_val(iter, N_dof, x_dummy, 0, val);

		bool mask[2] = {nb_bcond_iter_get_mask(iter, 0),
				nb_bcond_iter_get_mask(iter, 1)};
		val[0] *= factor;
		val[1] *= factor;
		nb_bcond_push(numeric_bcond, NB_NEUMANN, NB_BC_ON_POINT,
			      elem_id, mask, val);
	}
}

static void set_neumann_vtx(const nb_mesh2D_t *mesh,
			    uint32_t **elem_adj,
			    const nb_bcond_t *const bcond, 
			    double factor,
			    nb_bcond_t *numeric_bcond)
{
	uint8_t N_dof = nb_bcond_get_N_dof(bcond);
	uint16_t size = nb_bcond_iter_get_memsize();

	uint32_t memsize = size + N_dof * sizeof(double);
	char *memblock = nb_soft_allocate_mem(memsize);

	nb_bcond_iter_t *iter = (void*) memblock;
	double *val = (void*) (memblock + size);

	nb_bcond_iter_init(iter);
	nb_bcond_iter_set_conditions(iter, bcond, NB_NEUMANN,
				     NB_BC_ON_POINT);

	while (nb_bcond_iter_has_more(iter)) {
		nb_bcond_iter_go_next(iter);

		uint32_t model_node_id = nb_bcond_iter_get_id(iter);
		uint32_t node_id = nb_mesh2D_get_invtx(mesh, model_node_id);
		uint32_t subsgm_id[4];
		get_subsgm_adj_to_node(mesh, node_id, subsgm_id);

		double x[2];
		x[0] = nb_mesh2D_node_get_x(mesh, node_id);
		x[1] = nb_mesh2D_node_get_y(mesh, node_id);
		nb_bcond_iter_get_val(iter, N_dof, x, 0, val);

		bool mask[2] = {nb_bcond_iter_get_mask(iter, 0),
				nb_bcond_iter_get_mask(iter, 1)};

		set_neumann_subsgm_adj_to_node(mesh, elem_adj, N_dof,
					       factor, val, mask,
					       subsgm_id[0], subsgm_id[1],
					       numeric_bcond);

		set_neumann_subsgm_adj_to_node(mesh, elem_adj, N_dof,
					       factor, val, mask,
					       subsgm_id[2], subsgm_id[3],
					       numeric_bcond);
	}
	nb_bcond_iter_finish(iter);

	nb_soft_free_mem(memsize, memblock);
}

static void set_neumann_subsgm_adj_to_node(const nb_mesh2D_t *mesh,
					   uint32_t **elem_adj,
					   uint8_t N_dof,
					   double factor, double val[2],
					   bool mask[2], uint32_t sgm_id,
					   uint32_t subsgm_id,
					   nb_bcond_t *numeric_bcond)
{
	double subsgm_length =
		nb_mesh2D_insgm_subsgm_get_length(mesh, sgm_id, subsgm_id);
	double integral_factor = 0.5 * subsgm_length;
	uint32_t elem_id = elem_adj[sgm_id][subsgm_id];

	val[0] *= integral_factor * factor;
	val[1] *= integral_factor * factor;
	nb_bcond_push(numeric_bcond, NB_NEUMANN, NB_BC_ON_POINT,
		      elem_id, mask, val);
}

static void get_subsgm_adj_to_node(const nb_mesh2D_t *mesh, uint32_t node_id,
				   uint32_t subsgm_data[4])
{
	uint32_t N_sgm = nb_mesh2D_get_N_insgm(mesh);
	uint8_t cnt = 0;
	for (uint32_t i = 0; i < N_sgm; i++) {
		uint32_t N = nb_mesh2D_insgm_get_N_nodes(mesh, i);
		for (uint16_t j = 0; j < N; j++) {
			uint32_t nid = nb_mesh2D_insgm_get_node(mesh, i, j);
			if (nid == node_id) {
				if (j > 0) {
					subsgm_data[cnt * 2] = i;
					subsgm_data[cnt*2+1] = j - 1;
					cnt ++;
				}

				if (j < N - 1) {
					subsgm_data[cnt * 2] = i;
					subsgm_data[cnt*2+1] = j;
					cnt ++;
				}
				
				if (cnt < 2)
					break;
				else
					goto EXIT;
			}
		}
	}
EXIT:
	return;
}

static void set_dirichlet_sgm(const nb_mesh2D_t *mesh,
			      uint32_t **elem_adj,
			      const nb_bcond_t *const bcond, 
			      double factor,
			      nb_bcond_t *numeric_bcond)
{
	uint8_t N_dof = nb_bcond_get_N_dof(bcond);
	uint16_t size = nb_bcond_iter_get_memsize();
	nb_bcond_iter_t *iter = nb_allocate_on_stack(size);
	nb_bcond_iter_init(iter);
	nb_bcond_iter_set_conditions(iter, bcond, NB_DIRICHLET,
				     NB_BC_ON_SEGMENT);
	double *val = nb_allocate_on_stack(N_dof * sizeof(double));
	while (nb_bcond_iter_has_more(iter)) {
		nb_bcond_iter_go_next(iter);
		uint32_t sgm_id = nb_bcond_iter_get_id(iter);
		uint32_t N = nb_mesh2D_insgm_get_N_subsgm(mesh, sgm_id);
		for (uint32_t i = 0; i < N; i++) {
			uint32_t elem_id = elem_adj[sgm_id][i];
			
			bool mask[2] = {nb_bcond_iter_get_mask(iter, 0),
					nb_bcond_iter_get_mask(iter, 1)};
			
			double x[2] = {nb_mesh2D_elem_get_x(mesh, elem_id),
				       nb_mesh2D_elem_get_y(mesh, elem_id)};
			nb_bcond_iter_get_val(iter, N_dof, x, 0, val);
			val[0] *= factor;
			val[1] *= factor;

			nb_bcond_push(numeric_bcond, NB_DIRICHLET,
				      NB_BC_ON_POINT, elem_id, mask, val);
		}
	}
	nb_bcond_iter_finish(iter);	
}

static void set_dirichlet_vtx(const nb_mesh2D_t *mesh,
			      uint32_t **elem_adj,
			      const nb_bcond_t *const bcond, 
			      double factor,
			      nb_bcond_t *numeric_bcond)
{
	uint8_t N_dof = nb_bcond_get_N_dof(bcond);
	uint16_t size = nb_bcond_iter_get_memsize();
	nb_bcond_iter_t *iter = nb_allocate_on_stack(size);
	nb_bcond_iter_init(iter);
	nb_bcond_iter_set_conditions(iter, bcond, NB_DIRICHLET,
				     NB_BC_ON_POINT);
	double *val = nb_allocate_on_stack(N_dof * sizeof(*val));
	while (nb_bcond_iter_has_more(iter)) {
		nb_bcond_iter_go_next(iter);
		uint32_t model_node_id = nb_bcond_iter_get_id(iter);
		uint32_t node_id = nb_mesh2D_get_invtx(mesh, model_node_id);
		uint32_t subsgm_data[4];
		get_subsgm_adj_to_node(mesh, node_id, subsgm_data);
			
		bool mask[2] = {nb_bcond_iter_get_mask(iter, 0),
				nb_bcond_iter_get_mask(iter, 1)};
			
		double x[2] = {nb_mesh2D_node_get_x(mesh, node_id),
			       nb_mesh2D_node_get_y(mesh, node_id)};
		nb_bcond_iter_get_val(iter, N_dof, x, 0, val);
		val[0] *= factor;
		val[1] *= factor;

		uint8_t N_subsgm_adj_to_node = 2;
		for (uint8_t i = 0; i < N_subsgm_adj_to_node; i++) {
			uint32_t sgm_id = subsgm_data[i*2];
			uint32_t subsgm_id = subsgm_data[i*2+1];
			uint32_t elem_id = elem_adj[sgm_id][subsgm_id];
			nb_bcond_push(numeric_bcond, NB_DIRICHLET,
				      NB_BC_ON_POINT, elem_id, mask, val);
		}
	}
	nb_bcond_iter_finish(iter);
}

void nb_cvfa_set_numeric_bconditions(nb_sparse_t *K, double *F,
				     const nb_mesh2D_t *const mesh,
				     nb_bcond_t *bcond)
{
	uint8_t N_dof = nb_bcond_get_N_dof(bcond);
	uint16_t size = nb_bcond_iter_get_memsize();
	nb_bcond_iter_t *iter = nb_allocate_on_stack(size);

	nb_bcond_iter_init(iter);
	nb_bcond_iter_set_conditions(iter, bcond,
				     NB_NEUMANN,
				     NB_BC_ON_POINT);
	while (nb_bcond_iter_has_more(iter)) {
		nb_bcond_iter_go_next(iter);
		uint32_t elem_id = nb_bcond_iter_get_id(iter);
		set_numeric_bcond_neumann(mesh, elem_id, N_dof,
					  F, iter);
	}
	nb_bcond_iter_finish(iter);

	nb_bcond_iter_init(iter);
	nb_bcond_iter_set_conditions(iter, bcond,
				     NB_DIRICHLET,
				     NB_BC_ON_POINT);
	while (nb_bcond_iter_has_more(iter)) {
		nb_bcond_iter_go_next(iter);
		uint32_t elem_id = nb_bcond_iter_get_id(iter);
		set_numeric_bcond_dirichlet(mesh, K, N_dof,
					    F, iter, elem_id);
	}
	nb_bcond_iter_finish(iter);
}

static void set_numeric_bcond_neumann(const nb_mesh2D_t *mesh,
				      uint32_t elem_id, 
				      uint8_t N_dof, double* F,
				      nb_bcond_iter_t *iter)
{
	double x[2];
	x[0] = nb_mesh2D_elem_get_x(mesh, elem_id);
	x[1] = nb_mesh2D_elem_get_y(mesh, elem_id);

	double *val = nb_allocate_on_stack(N_dof * sizeof(*val));
	nb_bcond_iter_get_val(iter, N_dof, x, 0, val);
	for (uint8_t j = 0; j < N_dof; j++) {
		bool mask = nb_bcond_iter_get_mask(iter, j);
		if (mask) {
			uint32_t mtx_id = elem_id * N_dof + j;
			F[mtx_id] += val[j];
		}
	}
}

static void set_numeric_bcond_dirichlet(const nb_mesh2D_t *mesh,
					nb_sparse_t* K, uint8_t N_dof,
					double* F,
					nb_bcond_iter_t *iter,
					uint32_t elem_id)
{
	double x[2];
	x[0] = nb_mesh2D_elem_get_x(mesh, elem_id);
	x[1] = nb_mesh2D_elem_get_y(mesh, elem_id);

	double *val = nb_allocate_on_stack(N_dof * sizeof(*val));
	nb_bcond_iter_get_val(iter, N_dof, x, 0, val);
	for (uint8_t j = 0; j < N_dof; j++) {
		bool mask = nb_bcond_iter_get_mask(iter, j);
		if (mask) {
			uint32_t mtx_id = elem_id * N_dof + j;
			nb_sparse_set_Dirichlet_condition(K, F, mtx_id,
							   val[j]);
		}
	}
}

void nb_cvfa_get_elem_adj_to_model_node(const nb_mesh2D_t *mesh,
					uint32_t vtx_id,
					uint32_t elem_id[2])
{
	uint32_t **elem_adj = allocate_elem_adj(mesh);
	get_elem_adj(mesh, elem_adj);

	uint32_t node_id = nb_mesh2D_get_invtx(mesh, vtx_id);
	uint32_t subsgm_data[4];
	get_subsgm_adj_to_node(mesh, node_id, subsgm_data);
	
	uint32_t sgm_id = subsgm_data[0];
	uint32_t subsgm_id = subsgm_data[1];
	elem_id[0] = elem_adj[sgm_id][subsgm_id];
	
	sgm_id = subsgm_data[2];
	subsgm_id = subsgm_data[3];
	elem_id[1] = elem_adj[sgm_id][subsgm_id];

	free_elem_adj(elem_adj);
}
