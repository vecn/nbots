#ifndef __NB_PDE_BOT_FINITE_ELEMENT_SOLID_MECHANICS_STATIC_ELASTICITY2D_H__
#define __NB_PDE_BOT_FINITE_ELEMENT_SOLID_MECHANICS_STATIC_ELASTICITY2D_H__

#include <stdint.h>
#include <stdbool.h>
#include "nb/geometric_bot.h"

#include "nb/pde_bot/material.h"
#include "nb/pde_bot/common_solid_mechanics/analysis2D.h"
#include "nb/pde_bot/boundary_conditions/bcond.h"
#include "nb/pde_bot/finite_element/element.h"

int vcn_fem_compute_2D_Solid_Mechanics
			(const vcn_msh3trg_t *const mesh,
			 const vcn_fem_elem_t *const elemtype,
			 const nb_material_t *const material,
			 const nb_bcond_t *const bcond,
			 bool enable_self_weight,
			 double gravity[2],
			 nb_analysis2D_t analysis2D,
			 nb_analysis2D_params *params2D,
			 const bool* elements_enabled, /* NULL to enable all */
			 double* displacement, /* Output */
			 double* strain       /* Output */);

void vcn_fem_compute_stress_from_strain
			(uint32_t N_elements,
			 const vcn_fem_elem_t *const elemtype,
			 const nb_material_t *const material,
			 nb_analysis2D_t analysis2D,
			 double* strain,
			 const bool* elements_enabled /* NULL to enable all */,
			 double* stress /* Output */);

void vcn_fem_compute_von_mises(uint32_t N,
			       double *stress,
			       double *von_mises /* Output */);

#endif
