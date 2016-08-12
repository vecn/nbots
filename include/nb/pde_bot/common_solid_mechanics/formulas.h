#ifndef __NB_PDE_BOT_SOLID_MECHANICS_FORMULAS_H__
#define __NB_PDE_BOT_SOLID_MECHANICS_FORMULAS_H__

#include "nb/pde_bot/material.h"
#include "nb/pde_bot/common_solid_mechanics/analysis2D.h"

void nb_pde_get_constitutive_matrix(double D[4], 
				    const nb_material_t *material,
				    nb_analysis2D_t analysis2D);
double nb_pde_get_vm_stress(double sxx, double syy, double sxy);

#endif
