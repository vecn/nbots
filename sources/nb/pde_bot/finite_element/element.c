#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "nb/memory_bot.h"
#include "nb/pde_bot/finite_element/element.h"

#include "element_struct.h"

#define INV_3 0.33333333333333333333333333333

static nb_fem_elem_t *elem_allocate_mem(uint8_t N_nodes, uint8_t N_gp);
static void init_trg_linear(nb_fem_elem_t *elem);
static void init_quad_linear(nb_fem_elem_t *elem);

nb_fem_elem_t *nb_fem_elem_create(nb_elem_id type)
{
	nb_fem_elem_t *elem;
	switch(type) {
	case NB_TRG_LINEAR:
		elem = elem_allocate_mem(3, 1);
		init_trg_linear(elem);
		break;
	case NB_QUAD_LINEAR:
		elem = elem_allocate_mem(4, 4);
		init_quad_linear(elem);
		break;
	default:
		elem = elem_allocate_mem(3, 1);
		init_trg_linear(elem);
	}
	return elem;
}

static nb_fem_elem_t *elem_allocate_mem(uint8_t N_nodes, uint8_t N_gp)
{
	uint16_t w_memsize = N_gp * sizeof(double);
	uint16_t Ni_memsize = N_nodes * N_gp * sizeof(double);
	uint16_t elem_memsize = sizeof(nb_fem_elem_t);
	uint32_t memsize = elem_memsize + w_memsize + 3 * Ni_memsize;
	
	char *memblock = nb_allocate_mem(memsize);

	nb_fem_elem_t *elem = (void*) memblock;
	elem->gp_weight = (void*) (memblock + elem_memsize);
	elem->Ni = (void*) (memblock + elem_memsize + w_memsize);
	elem->dNi_dpsi = (void*) (memblock + elem_memsize +
				  w_memsize + Ni_memsize);
	elem->dNi_deta = (void*) (memblock + elem_memsize +
				  w_memsize + 2 * Ni_memsize);
	return elem;
}

static void init_trg_linear(nb_fem_elem_t *elem)
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

static void init_quad_linear(nb_fem_elem_t *elem)
{
	elem->type = NB_QUAD_LINEAR;
	elem->N_nodes = 4;
	elem->N_gp = 4;

	elem->gp_weight[0] = 1.0;
	elem->gp_weight[1] = 1.0;
	elem->gp_weight[2] = 1.0;
	elem->gp_weight[3] = 1.0;

	double Ni[16] = 
		{0.622008467928, 0.166666666667,
		 0.044658198739, 0.166666666667,
		 0.166666666667, 0.622008467928,
		 0.166666666667, 0.044658198739,
		 0.044658198739, 0.166666666667,
		 0.622008467928, 0.166666666667,
		 0.166666666667, 0.044658198739,
		 0.166666666667, 0.622008467928};

	memcpy(elem->Ni, Ni, 16 * sizeof(*(elem->Ni)));

	double dNi_dpsi[16] =
		{-0.394337567297, -0.394337567297,
		 -0.105662432703, -0.105662432703,
		 0.394337567297, 0.394337567297,
		 0.105662432703, 0.105662432703,
		 0.105662432703, 0.105662432703,
		 0.394337567297, 0.394337567297,
		 -0.105662432703, -0.105662432703,
		 -0.394337567297, -0.394337567297};

	memcpy(elem->dNi_dpsi, dNi_dpsi, 16 * sizeof(*(elem->dNi_dpsi)));

	double dNi_deta[16] =
		{-0.394337567297, -0.105662432703,
		 -0.105662432703, -0.394337567297,
		 -0.105662432703, -0.394337567297,
		 -0.394337567297, -0.105662432703,
		 0.105662432703, 0.394337567297,
		 0.394337567297, 0.105662432703,
		 0.394337567297, 0.105662432703,
		 0.105662432703, 0.394337567297};

	memcpy(elem->dNi_deta, dNi_deta, 16 * sizeof(*(elem->dNi_deta)));
}

void nb_fem_elem_destroy(nb_fem_elem_t* elem)
{
	nb_free_mem(elem);
}

uint8_t nb_fem_elem_get_N_gpoints(const nb_fem_elem_t *const elem)
{
	return elem->N_gp;
}

uint8_t nb_fem_elem_get_N_nodes(const nb_fem_elem_t *const elem)
{
	return elem->N_nodes;
}

double nb_fem_elem_weight_gp(const nb_fem_elem_t *const elem,
			      uint8_t gp_id)
{
	return elem->gp_weight[gp_id];	
}

double nb_fem_elem_Ni(const nb_fem_elem_t *const elem,
		       uint8_t node_id, uint8_t gp_id)
{
	return elem->Ni[node_id * elem->N_gp + gp_id];
}

double nb_fem_elem_dNi_dpsi(const nb_fem_elem_t *const elem,
			     uint8_t node_id, uint8_t gp_id)
{
	return elem->dNi_dpsi[node_id * elem->N_gp + gp_id];
}

double nb_fem_elem_dNi_deta(const nb_fem_elem_t *const elem,
			   uint8_t node_id, uint8_t gp_id)
{
	return elem->dNi_deta[node_id * elem->N_gp + gp_id];
}
