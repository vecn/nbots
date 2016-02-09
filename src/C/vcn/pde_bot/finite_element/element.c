#include <stdio.h>
#include <stdlib.h>

#include "vcn/pde_bot/finite_element/element.h"
#include "vcn/pde_bot/finite_element/elements/trg_linear.h"

#include "element_struct.h"

vcn_fem_elem_t* vcn_fem_elem_create(vcn_elem_id id)
{
	vcn_fem_elem_t *elem = calloc(1, sizeof(*elem));
	switch(id) {
	case VCN_TRG_LINEAR:
		vcn_fem_elem_init_trg_linear(elem);
	default:
		vcn_fem_elem_init_trg_linear(elem);
	}
	return elem;
}

void vcn_fem_elem_destroy(vcn_fem_elem_t* elemtype)
{
	free(elemtype->Ni);
	free(elemtype->dNi_dpsi);
	free(elemtype->dNi_deta);
	free(elemtype->psi);
	free(elemtype->eta);
	free(elemtype->gp_weight);
	free(elemtype);
}

inline uint32_t vcn_fem_elem_get_N_nodes(const vcn_fem_elem_t *const elemtype)
{
	return elemtype->N_nodes;
}

inline double vcn_fem_elem_eval_shape_function
		(const vcn_fem_elem_t *const elemtype, uint32_t node_id,
		 double psi, double eta)
{
	return elemtype->Ni[node_id](psi, eta);
}

inline uint32_t vcn_fem_elem_get_closest_Gauss_Point_to_the_ith_node
		(const vcn_fem_elem_t *const elemtype, uint32_t i)
{
	if (3 == elemtype->N_nodes && 
	    1 == elemtype->N_Gauss_points){
		/* Triangle of linear interpolation
		 *              o
		 *             / \
		 *            / * \    <---- Gauss Point
		 *           /_____\
		 *          o       o  <---- Nodes
		 */
		return 0;
	}
	printf("FEM Error: Unknown element type with %i nodes and %i Gauss points.\n",
	       elemtype->N_nodes, elemtype->N_Gauss_points);
	return 0;
}
