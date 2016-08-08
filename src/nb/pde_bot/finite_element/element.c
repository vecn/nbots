#include <stdio.h>
#include <stdlib.h>

#include "nb/pde_bot/finite_element/element.h"

#include "element_struct.h"
#include "elements/trg_linear.h"

vcn_fem_elem_t* vcn_fem_elem_create(vcn_elem_id type)
{
	vcn_fem_elem_t *elem = calloc(1, sizeof(*elem));
	elem->type = type;
	switch(type) {
	case NB_TRG_LINEAR:
		trg_linear_init(elem);
		break;
	default:
		trg_linear_init(elem);
		elem->type = NB_TRG_LINEAR;
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

uint8_t vcn_fem_elem_get_N_gpoints(const vcn_fem_elem_t *const elemtype)
{
	return elemtype->N_Gauss_points;
}

uint8_t vcn_fem_elem_get_N_nodes(const vcn_fem_elem_t *const elemtype)
{
	return elemtype->N_nodes;
}

double vcn_fem_elem_eval_shape_function
		(const vcn_fem_elem_t *const elemtype, uint8_t node_id,
		 double psi, double eta)
{
	return elemtype->Ni[node_id](psi, eta);
}

double vcn_fem_elem_eval_shape_function_on_gp
		(const vcn_fem_elem_t *const elemtype,
		 uint8_t node_id, uint8_t gp_id)
{
	double psi = elemtype->psi[gp_id];
	double eta = elemtype->eta[gp_id];
	return elemtype->Ni[node_id](psi, eta);
}
