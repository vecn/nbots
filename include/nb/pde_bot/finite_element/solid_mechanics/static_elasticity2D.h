#ifndef __NB_PDE_BOT_FINITE_ELEMENT_SOLID_MECHANICS_STATIC_ELASTICITY2D_H__
#define __NB_PDE_BOT_FINITE_ELEMENT_SOLID_MECHANICS_STATIC_ELASTICITY2D_H__

#include <stdint.h>
#include <stdbool.h>
#include "nb/eigen_bot.h"
#include "nb/geometric_bot.h"

#include "nb/pde_bot/material.h"
#include "nb/pde_bot/boundary_conditions.h"
#include "nb/pde_bot/finite_element/element.h"

int vcn_fem_compute_2D_Solid_Mechanics
			(const vcn_msh3trg_t *const mesh,
			 const vcn_fem_elem_t *const elemtype,
			 const vcn_fem_material_t *const material,
			 const vcn_bcond_t *const bmeshcond,
			 char enable_self_weight,
			 double gravity[2],
			 bool enable_plane_stress,
			 double thickness,
			 const bool* elements_enabled, /* NULL to enable all */
			 double* displacement, /* Output */
			 double* strain,       /* Output */
			 const char* logfile /* NULL if not required */);

void vcn_fem_compute_stress_from_strain
			(uint32_t N_elements,
			 uint32_t* elements_connectivity_matrix, 
			 const vcn_fem_elem_t *const elemtype,
			 const vcn_fem_material_t *const material,
			 bool enable_plane_stress,
			 double* strain,
			 const bool* elements_enabled /* NULL to enable all */,
			 double* stress /* Output */);

#endif
