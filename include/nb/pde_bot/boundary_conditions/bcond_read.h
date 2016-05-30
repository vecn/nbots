#ifndef __NB_PDE_BOT_BOUNDARY_CONDITIONS_BCOND_READ_H__
#define __NB_PDE_BOT_BOUNDARY_CONDITIONS_BCOND_READ_H__

#include <stdint.h>

#include "nb/pde_bot/boundary_conditions/bcond.h"

#include "nb/cfreader_cat.h"

int nb_bcond_read(nb_bcond_t *bcond, vcn_cfreader_t *cfreader);

#endif
