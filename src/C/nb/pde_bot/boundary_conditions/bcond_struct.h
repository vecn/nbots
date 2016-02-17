#ifndef __NB_PDE_BOT_BOUNDARY_CONDITIONS_BCOND_STRUCT_H__
#define __NB_PDE_BOT_BOUNDARY_CONDITIONS_BCOND_STRUCT_H__

#include <stdint.h>

#include "nb/container_bot.h"

struct nb_bcond_s {
	uint8_t N_dof; /* Degrees of freedom */
	nb_container_t *dirichlet_vtx;
	nb_container_t *neumann_vtx;
	nb_container_t *dirichlet_sgm;
	nb_container_t *neumann_sgm;
};

#endif
