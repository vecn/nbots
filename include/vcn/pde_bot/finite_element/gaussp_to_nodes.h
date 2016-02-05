#ifndef __VCN_PDE_BOT_FINITE_ELEMENT_GAUSSP_TO_NODES_H__
#define __VCN_PDE_BOT_FINITE_ELEMENT_GAUSSP_TO_NODES_H__

#include <stdint.h>
#include "vcn/geometric_bot/mesh/elements/triangle.h"
#include "vcn/pde_bot/finite_element/element.h"

void vcn_fem_interpolate_from_Gauss_points_to_nodes
			(const vcn_msh3trg_t *const mesh,
			 const vcn_fem_elem_t *const elemtype,
			 uint32_t N_components,
			 double* values_on_GP_from_elements,
			 double* values_interpolated_on_nodes /* Output */);

#endif
