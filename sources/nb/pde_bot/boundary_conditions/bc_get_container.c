#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "nb/container_bot.h"
#include "nb/pde_bot/boundary_conditions/bcond.h"

#include "bcond_struct.h"
#include "bc_get_container.h"

static nb_container_t* get_dirichlet_container(const nb_bcond_t *bcond,
					       nb_bcond_where type_elem);
static nb_container_t* get_neumann_container(const nb_bcond_t *bcond,
					     nb_bcond_where type_elem);



nb_container_t *nb_bcond_get_container(const nb_bcond_t *const bcond,
					nb_bcond_id type_id,
					nb_bcond_where type_elem)
{
	nb_container_t *container;
	switch (type_id) {
	case NB_DIRICHLET:
		container = get_dirichlet_container(bcond, type_elem);
		break;
	case NB_NEUMANN:
		container = get_neumann_container(bcond, type_elem);
		break;
	default:
		container = NULL;
	}
	return container;
}

static nb_container_t* get_dirichlet_container(const nb_bcond_t *bcond,
					       nb_bcond_where type_elem)
{
	nb_container_t *container;
	switch (type_elem) {
	case NB_BC_ON_POINT:
		container = bcond->dirichlet_vtx;
		break;
	case NB_BC_ON_SEGMENT:
		container = bcond->dirichlet_sgm;
		break;
	default:
		container = NULL;
	}
	return container;
}

static nb_container_t* get_neumann_container(const nb_bcond_t *bcond,
					     nb_bcond_where type_elem)
{
	nb_container_t *container;
	switch (type_elem) {
	case NB_BC_ON_POINT:
		container = bcond->neumann_vtx;
		break;
	case NB_BC_ON_SEGMENT:
		container = bcond->neumann_sgm;
		break;
	default:
		container = NULL;
	}
	return container;
}
