#ifndef __NB_PDE_BOT_FINITE_ELEMENT_GAUSSP_TO_NODES_H__
#define __NB_PDE_BOT_FINITE_ELEMENT_GAUSSP_TO_NODES_H__

#include <stdint.h>

#include "nb/geometric_bot.h"
#include "nb/pde_bot/finite_element/element.h"

int nb_fem_interpolate_from_gpoints_to_nodes
			(const nb_partition_t *const part,
			 const nb_fem_elem_t *const elem,
			 uint32_t N_comp,
			 const double* gp_values,
			 double* nodal_values /* Output */);

void nb_fem_get_error_on_gpoints(const nb_partition_t *const part,
				  const nb_fem_elem_t *const elem,
				  uint32_t N_comp,
				  double* gp_values,
				  double* gp_error);
#endif
