#ifndef __NB_PDE_BOT_FINITE_ELEMENT_SOLID_MECHANICS_STATIC_DAMAGE2D_H__
#define __NB_PDE_BOT_FINITE_ELEMENT_SOLID_MECHANICS_STATIC_DAMAGE2D_H__

#include <stdint.h>
#include <stdbool.h>
#include "nb/geometric_bot.h"

#include "nb/pde_bot/material.h"
#include "nb/pde_bot/common_solid_mechanics/analysis2D.h"
#include "nb/pde_bot/boundary_conditions/bcond.h"
#include "nb/pde_bot/finite_element/element.h"

typedef struct nb_fem_implicit_s nb_fem_implicit_t;

nb_fem_implicit_t* nb_fem_implicit_create();
void nb_fem_implicit_destroy(nb_fem_implicit_t* isparams);

void nb_fem_implicit_set_N_steps(nb_fem_implicit_t* isparams,
				  uint32_t N_time_steps);
void nb_fem_implicit_set_N_max_iter(nb_fem_implicit_t* isparams,
				     uint32_t N_max_iter);
void nb_fem_implicit_set_N_max_iter_without_enhance
(nb_fem_implicit_t* isparams, uint32_t N_max_iter);
void nb_fem_implicit_set_residual_tolerance
(nb_fem_implicit_t* isparams, double penergy_tol);

uint32_t nb_fem_implicit_get_N_steps(nb_fem_implicit_t* isparams);
uint32_t nb_fem_implicit_get_N_max_iter(nb_fem_implicit_t* isparams);
uint32_t nb_fem_implicit_get_N_max_iter_without_enhance
(nb_fem_implicit_t* isparams);
double nb_fem_implicit_get_residual_tolerance
(nb_fem_implicit_t* isparams);

void nb_fem_compute_2D_Non_Linear_Solid_Mechanics
		(const nb_partition_t *const part,
		 const nb_fem_elem_t *const elemtype,
		 const nb_material_t *const material,
		 const nb_bcond_t *const bcond,
		 bool enable_self_weight,
		 double gravity[2],
		 bool enable_Cholesky_solver,
		 nb_analysis2D_t analysis2D,
		 nb_analysis2D_params *params2D,
		 nb_fem_implicit_t* params,
		 const char* logfile);


#endif
