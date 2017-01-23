#ifndef __NB_PDE_BOT_CONTROL_VOLUME_SOLID_MECHANICS_BCONDITIONS_H__
#define __NB_PDE_BOT_CONTROL_VOLUME_SOLID_MECHANICS_BCONDITIONS_H__

#include "nb/solver_bot.h"
#include "nb/geometric_bot.h"
#include "nb/pde_bot/boundary_conditions/bcond.h"

void nb_cvfa_set_bconditions(const nb_mesh2D_t *mesh,
			     const nb_material_t *material,
			     nb_analysis2D_t analysis2D,
			     nb_sparse_t* K, double* F, 
			     const nb_bcond_t *const bcond,
			     double factor);

void nb_cvfa_get_numeric_bconditions(const nb_mesh2D_t *mesh,
				     const nb_material_t *material,
				     nb_analysis2D_t analysis2D,
				     const nb_bcond_t *const bcond,
				     double factor,
				     nb_bcond_t *numeric_bcond);

void nb_cvfa_set_numeric_bconditions(nb_sparse_t *K, double *F,
				     const nb_mesh2D_t *const mesh,
				     const nb_bcond_t *bcond);

void nb_cvfa_get_elem_adj_to_model_node(const nb_mesh2D_t *mesh,
					uint32_t vtx_id,
					uint32_t elem_id[2]);
#endif
