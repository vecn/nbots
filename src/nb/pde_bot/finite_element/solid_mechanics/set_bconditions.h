#ifndef __NB_PDE_BOT_FINITE_ELEMENT_SOLID_MECHANICS_BCONDITIONS_H__
#define __NB_PDE_BOT_FINITE_ELEMENT_SOLID_MECHANICS_BCONDITIONS_H__

#include "nb/eigen_bot.h"
#include "nb/geometric_bot.h"
#include "nb/pde_bot/boundary_conditions/bcond.h"
#include "nb/pde_bot/boundary_conditions/bcond_iter.h"

void nb_set_bconditions(const nb_partition_t *part,
			vcn_sparse_t* K, double* F, 
			const nb_bcond_t *const bcond,
			double factor);

#endif
