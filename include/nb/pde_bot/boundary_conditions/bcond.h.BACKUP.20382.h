#ifndef __NB_PDE_BOT_BOUNDARY_CONDITIONS_H__
#define __NB_PDE_BOT_BOUNDARY_CONDITIONS_H__

#include <stdint.h>

<<<<<<< HEAD
#include "vcn/container_bot.h"
=======
#include "vcn/cfreader_cat.h"
>>>>>>> 0a4b4219cfcc0385e64339b520a04ee71c142a2a

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
<<<<<<< HEAD
void nb_bcond_init(void *bcond_ptr);
=======
void nb_bcond_init(void *bcond_ptr, uint8_t N_dof);
>>>>>>> 0a4b4219cfcc0385e64339b520a04ee71c142a2a
void nb_bcond_copy(void *bcond_ptr, const void *const src_bcond_ptr);
void nb_bcond_clear(void *bcond_ptr);
void* nb_bcond_create(uint8_t N_dof);
void* nb_bcond_clone(const void *const bcond_ptr);
void nb_bcond_destroy(void *bcond_ptr);

uint8_t nb_bcond_get_N_dof(const nb_bcond_t *const bcond);
<<<<<<< HEAD
void nb_bcond_read(nb_bcond_t *bcond, vcn_cfreader_t *cfreader);
=======
int nb_bcond_read(nb_bcond_t *bcond, vcn_cfreader_t *cfreader);
>>>>>>> 0a4b4219cfcc0385e64339b520a04ee71c142a2a
void nb_bcond_push(nb_bcond_t *bcond, nb_bcond_id type_id,
		   nb_bcond_where type_elem, uint32_t elem_id,
		   const bool dof_mask[], const double value[]);

#endif
