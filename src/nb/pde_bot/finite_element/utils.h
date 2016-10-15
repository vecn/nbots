#ifndef __NB_PDE_BOT_FINITE_ELEMENT_UTILS_H__
#define __NB_PDE_BOT_FINITE_ELEMENT_UTILS_H__

#include "nb/geometric_bot.h"
#include "nb/pde_bot/material.h"
#include "nb/pde_bot/finite_element/element.h"

void nb_fem_get_derivatives(const nb_fem_elem_t *elem,
			    int gp_id, double Jinv[4],
			    double *dNi_dx, double *dNi_dy);

bool nb_fem_elem_is_distorted(double detJ);

double nb_fem_get_jacobian(const nb_fem_elem_t *elem, uint32_t id,
			   const nb_mesh2D_t *part, int gp_id,
			   double Jinv[4]);

#endif
