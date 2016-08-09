#include <stdio.h>
#include <stdlib.h>

#include "nb/pde_bot/finite_element/element.h"

#include "element_struct.h"

#define INV_3 0.33333333333333333333333333333

static vcn_fem_elem_t *elem_malloc(uint8_t N_nodes, uint8_t N_gp);
static void init_trg_linear(vcn_fem_elem_t *elem);

vcn_fem_elem_t *vcn_fem_elem_create(vcn_elem_id type)
{
	vcn_fem_elem_t *elem;
	switch(type) {
	case NB_TRG_LINEAR:
		elem = elem_malloc(3, 1);
		init_trg_linear(elem);
		break;
	default:
		elem = elem_malloc(3, 1);
		init_trg_linear(elem);
	}
	return elem;
}

static vcn_fem_elem_t *elem_malloc(uint8_t N_nodes, uint8_t N_gp)
{
	uint16_t w_memsize = N_gp * sizeof(double);
	uint16_t Ni_memsize = N_nodes * N_gp * sizeof(double);
	uint16_t elem_memsize = sizeof(vcn_fem_elem_t);
	uint32_t memsize = elem_memsize + w_memsize + 3 * Ni_memsize;
	
	char *memblock = malloc(memsize);

	vcn_fem_elem_t *elem = (void*) memblock;
	elem->gp_weight = (void*) (memblock + elem_memsize);
	elem->Ni = (void*) (memblock + elem_memsize + w_memsize);
	elem->dNi_dpsi = (void*) (memblock + elem_memsize +
				  w_memsize + Ni_memsize);
	elem->dNi_deta = (void*) (memblock + elem_memsize +
				  w_memsize + 2 * Ni_memsize);
	return elem;
}

static void init_trg_linear(vcn_fem_elem_t *elem)
{
	elem->type = NB_TRG_LINEAR;
	elem->N_nodes = 3;
	elem->N_gp = 1;

	elem->gp_weight[0] = 0.5;

	elem->Ni[0] = INV_3;
	elem->Ni[1] = INV_3;
	elem->Ni[2] = INV_3;

	elem->dNi_dpsi[0] = -1.0;
	elem->dNi_dpsi[1] = 1.0;
	elem->dNi_dpsi[2] = 0.0;

	elem->dNi_deta[0] = -1.0;
	elem->dNi_deta[1] = 0.0;
	elem->dNi_deta[2] = 1.0;
}

void vcn_fem_elem_destroy(vcn_fem_elem_t* elem)
{
	free(elem);
}

uint8_t vcn_fem_elem_get_N_gpoints(const vcn_fem_elem_t *const elem)
{
	return elem->N_gp;
}

uint8_t vcn_fem_elem_get_N_nodes(const vcn_fem_elem_t *const elem)
{
	return elem->N_nodes;
}

double vcn_fem_elem_weight_gp(const vcn_fem_elem_t *const elem,
			      uint8_t gp_id)
{
	return elem->gp_weight[gp_id];	
}

double vcn_fem_elem_Ni(const vcn_fem_elem_t *const elem,
		       uint8_t node_id, uint8_t gp_id)
{
	return elem->Ni[node_id * elem->N_gp + gp_id];
}

double vcn_fem_elem_dNi_dpsi(const vcn_fem_elem_t *const elem,
			     uint8_t node_id, uint8_t gp_id)
{
	return elem->dNi_dpsi[node_id * elem->N_gp + gp_id];
}

double vcn_fem_elem_dNi_deta(const vcn_fem_elem_t *const elem,
			   uint8_t node_id, uint8_t gp_id)
{
	return elem->dNi_deta[node_id * elem->N_gp + gp_id];
}
