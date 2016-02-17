#ifndef __NB_PDE_BOT_BOUNDARY_CONDITIONS_H__
#define __NB_PDE_BOT_BOUNDARY_CONDITIONS_H__

#include <stdint.h>

#include "nb/cfreader_cat.h"

typedef enum {
	NB_DIRICHLET,
	NB_NEUMANN
} nb_bcond_id;

typedef enum {
	NB_BC_ON_POINT,
	NB_BC_ON_SEGMENT
} nb_bcond_where;

typedef struct nb_bcond_s nb_bcond_t;

uint16_t nb_bcond_get_memsize(uint8_t N_dof);
void nb_bcond_init(void *bcond_ptr, uint8_t N_dof);
void nb_bcond_copy(void *bcond_ptr, const void *const src_bcond_ptr);
void nb_bcond_finish(void *bcond_ptr);
void* nb_bcond_create(uint8_t N_dof);
void* nb_bcond_clone(const void *const bcond_ptr);
void nb_bcond_destroy(void *bcond_ptr);
void nb_bcond_clear(void *bcond_ptr);

uint8_t nb_bcond_get_N_dof(const nb_bcond_t *const bcond);
int nb_bcond_read(nb_bcond_t *bcond, vcn_cfreader_t *cfreader);
void nb_bcond_push(nb_bcond_t *bcond, nb_bcond_id type_id,
		   nb_bcond_where type_elem, uint32_t elem_id,
		   const bool dof_mask[], const double value[]);

#endif
