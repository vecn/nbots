#ifndef __NB_PDE_BOT_FINITE_ELEMENT_SOLID_MECHANICS_STATIC_DAMAGE2D_H__
#define __NB_PDE_BOT_FINITE_ELEMENT_SOLID_MECHANICS_STATIC_DAMAGE2D_H__

#include <stdint.h>
#include <stdbool.h>
#include "nb/geometric_bot.h"

#include "nb/pde_bot/material.h"
#include "nb/pde_bot/common_solid_mechanics/analysis2D.h"
#include "nb/pde_bot/boundary_conditions/bcond.h"
#include "nb/pde_bot/finite_element/element.h"

typedef struct vcn_fem_implicit_s vcn_fem_implicit_t;

vcn_fem_implicit_t* vcn_fem_implicit_create();
void vcn_fem_implicit_destroy(vcn_fem_implicit_t* isparams);

void vcn_fem_implicit_set_N_steps(vcn_fem_implicit_t* isparams,
				  uint32_t N_time_steps);
void vcn_fem_implicit_set_N_max_iter(vcn_fem_implicit_t* isparams,
				     uint32_t N_max_iter);
void vcn_fem_implicit_set_N_max_iter_without_enhance
(vcn_fem_implicit_t* isparams, uint32_t N_max_iter);
void vcn_fem_implicit_set_residual_tolerance
(vcn_fem_implicit_t* isparams, double penergy_tol);

uint32_t vcn_fem_implicit_get_N_steps(vcn_fem_implicit_t* isparams);
uint32_t vcn_fem_implicit_get_N_max_iter(vcn_fem_implicit_t* isparams);
uint32_t vcn_fem_implicit_get_N_max_iter_without_enhance
(vcn_fem_implicit_t* isparams);
double vcn_fem_implicit_get_residual_tolerance
(vcn_fem_implicit_t* isparams);

void vcn_fem_compute_2D_Non_Linear_Solid_Mechanics
(const vcn_msh3trg_t *const mesh,
 const vcn_fem_elem_t *const elemtype,
 const nb_material_t *const material,
 const nb_bcond_t *const bcond,
 bool enable_self_weight,
 double gravity[2],
 bool enable_Cholesky_solver,
 nb_analysis2D_t analysis2D,
 nb_analysis2D_params *params2D,
 vcn_fem_implicit_t* params,
 bool restore_computation, /* Restore computation after crash */
 const char* logfile);


#endif
