#include <stdio.h>
#include <stdlib.h>

#include "nb/pde_bot/finite_element/element.h"

#include "element_struct.h"

#define INV_3 (0.33333333333333333333333333333)

static vcn_fem_elem_t *trg_linear_malloc(void)
static void init_trg_linear(vcn_fem_elem_t *elem);

vcn_fem_elem_t *vcn_fem_elem_create(vcn_elem_id type)
{
	vcn_fem_elem_t *elem;
	elem->type = type;
	switch(type) {
	case NB_TRG_LINEAR:
		elem = trg_linear_malloc();
		trg_linear_init(elem);
		break;
	default:
		elem = trg_linear_malloc();
		trg_linear_init(elem);
		elem->type = NB_TRG_LINEAR;
	}
	return elem;
}

static vcn_fem_elem_t *trg_linear_malloc(void)
{
	/* AQUI VOY */
}

static void init_trg_linear(vcn_fem_elem_t *elem)
{
	

	elem->N_nodes = 3;
	elem->N_gp = 1;

	elem->Ni = calloc(elem->N_nodes, sizeof(*(elem->Ni)));
	elem->dNi_dpsi = calloc(elem->N_nodes,
				    sizeof(*(elem->dNi_dpsi)));
	elem->dNi_deta = calloc(elem->N_nodes,
				    sizeof(*(elem->dNi_deta)));  
	elem->Ni[0] = N1;
	elem->Ni[1] = N2;
	elem->Ni[2] = N3;
	elem->dNi_dpsi[0] = dN1_dpsi;
	elem->dNi_dpsi[1] = dN2_dpsi;
	elem->dNi_dpsi[2] = dN3_dpsi;
	elem->dNi_deta[0] = dN1_deta;
	elem->dNi_deta[1] = dN2_deta;
	elem->dNi_deta[2] = dN3_deta;
  
	elem->psi = malloc(elem->N_Gauss_points * sizeof(*(elem->psi)));
	elem->eta = malloc(elem->N_Gauss_points * sizeof(*(elem->eta)));
	elem->gp_weight = malloc(elem->N_Gauss_points *
				 sizeof(*(elem->gp_weight)));

	elem->psi[0] = INV_3;
	elem->eta[0] = INV_3;

	elem->gp_weight[0] = 0.5;
}

void vcn_fem_elem_destroy(vcn_fem_elem_t* elemtype)
{
	free(elemtype);
}

uint8_t vcn_fem_elem_get_N_gpoints(const vcn_fem_elem_t *const elemtype)
{
	return elemtype->N_gp;
}

uint8_t vcn_fem_elem_get_N_nodes(const vcn_fem_elem_t *const elemtype)
{
	return elemtype->N_nodes;
}

double vcn_fem_elem_Ni(const vcn_fem_elem_t *const elemtype,
		       uint8_t node_id, uint8_t gp_id)
{
	uint8_t N_gp = elemtype->N_Gauss_points;
	return elemtype->Ni[node_id * N_gp + gp_id];
}

double vcn_fem_elem_dNi_dpsi(const vcn_fem_elem_t *const elemtype,
			     uint8_t node_id, uint8_t gp_id)
{
	uint8_t N_gp = elemtype->N_Gauss_points;
	return elemtype->dNi_dpsi[node_id * N_gp + gp_id];
}

double vcn_fem_elem_dNi_deta(const vcn_fem_elem_t *const elemtype,
			   uint8_t node_id, uint8_t gp_id)
{
	uint8_t N_gp = elemtype->N_Gauss_points;
	return elemtype->dNi_deta[node_id * N_gp + gp_id];
}
