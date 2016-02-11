#ifndef __NB_PDE_BOT_BOUNDARY_CONDITIONS_H__
#define __NB_PDE_BOT_BOUNDARY_CONDITIONS_H__

#include <stdint.h>

#include "vcn/container_bot.h"

typedef enum {
	NB_DIRICHLET,
	NB_NEUMANN
} nb_bcond_id;

typedef enum {
	NB_BC_ON_POINT,
	NB_BC_ON_SEGMENT
} nb_bcond_where;

typedef struct nb_bcond_s nb_bcond_t;

nb_bcond_t* nb_bcond_create(uint8_t N_dof);
nb_bcond_t* nb_bcond_clone(const nb_bcond_t *const bcond);
nb_bcond_t* nb_bcond_read(vcn_cfreader_t *cfreader);
void nb_bcond_save(const nb_bcond_t *const bcond,
		   const char* filename);
void nb_bcond_copy(nb_bcond_t* dest, const nb_bcond_t *const src);
void nb_bcond_clear(nb_bcond_t* bcond);
void nb_bcond_destroy(nb_bcond_t* bcond);
void nb_bcond_push(nb_bcond_t *bcond, nb_bcond_id type_id,
		   nb_bcond_where type_elem, uint32_t elem_id,
		   bool dof_mask[], double value[]);

#endif
