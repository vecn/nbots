#ifndef __NB_PDE_BOT_CONTROL_VOLUME_SOLID_MECHANICS_STATIC_ELASTICITY2D_H__
#define __NB_PDE_BOT_CONTROL_VOLUME_SOLID_MECHANICS_STATIC_ELASTICITY2D_H__

#include <stdint.h>
#include <stdbool.h>
#include "nb/geometric_bot.h"

#include "nb/pde_bot/material.h"
#include "nb/pde_bot/common_solid_mechanics/analysis2D.h"
#include "nb/pde_bot/boundary_conditions/bcond.h"
#include "nb/pde_bot/control_volume/solid_mechanics/static_elasticity2D.h"

int nb_cvfa_compute_2D_Solid_Mechanics
			(const nb_partition_t *const part,
			 const nb_material_t *const material,
			 const nb_bcond_t *const bcond,
			 bool enable_self_weight,
			 double gravity[2],
			 nb_analysis2D_t analysis2D,
			 nb_analysis2D_params *params2D,
			 double *displacement, /* Output */
			 double *strain       /* Output */);

void nb_cvfa_compute_stress_from_strain(const nb_partition_t *part,
					const nb_material_t *const material,
					nb_analysis2D_t analysis2D,
					double* strain,
					double* stress /* Output */);

#endif
