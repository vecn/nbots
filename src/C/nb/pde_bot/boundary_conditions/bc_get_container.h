#ifndef __NB_PDE_BOT_BOUNDARY_CONDITIONS_GET_CONTAINER_H__
#define __NB_PDE_BOT_BOUNDARY_CONDITIONS_GET_CONTAINER_H__

#include "nb/container_bot.h"
#include "nb/pde_bot/boundary_conditions/bcond.h"

vcn_container_t *nb_bcond_get_container(const nb_bcond_t *const bcond,
					nb_bcond_id type_id,
					nb_bcond_where type_elem);

#endif
