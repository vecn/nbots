#ifndef __NB_PDE_BOT_FINITE_ELEMENT_GAUSSP_TO_NODES_H__
#define __NB_PDE_BOT_FINITE_ELEMENT_GAUSSP_TO_NODES_H__

#include <stdint.h>
#include "nb/geometric_bot/mesh/elements2D/triangles.h"
#include "nb/pde_bot/finite_element/element.h"

int vcn_fem_interpolate_from_gpoints_to_nodes
			(const vcn_msh3trg_t *const mesh,
			 const vcn_fem_elem_t *const elem,
			 uint32_t N_comp,
			 const double* gp_values,
			 double* nodal_values /* Output */);

void vcn_fem_get_error_on_gpoints(const vcn_msh3trg_t *const mesh,
				  const vcn_fem_elem_t *const elem,
				  uint32_t N_comp,
				  double* gp_values,
				  double* gp_error);
#endif
