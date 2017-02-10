#ifndef __NB_PDE_BOT_FINITE_ELEMENT_SOLID_MECHANICS_BCONDITIONS_H__
#define __NB_PDE_BOT_FINITE_ELEMENT_SOLID_MECHANICS_BCONDITIONS_H__

#include "nb/solver_bot.h"
#include "nb/geometric_bot.h"
#include "nb/pde_bot/boundary_conditions/bcond.h"

void nb_fem_set_bconditions(const nb_mesh2D_t *part,
			    nb_sparse_t* K, double* F, 
			    const nb_bcond_t *const bcond,
			    double factor);
void nb_fem_set_dirichlet_bconditions(const nb_mesh2D_t *part,
			    nb_sparse_t* K, double* F, 
			    const nb_bcond_t *const bcond,
			    double factor);

#endif
