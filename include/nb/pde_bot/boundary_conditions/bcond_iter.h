#ifndef __NB_PDE_BOT_BCONDITIONS_ITERATOR_H__
#define __NB_PDE_BOT_BCONDITIONS_ITERATOR_H__

#include <stdint.h>
#include <stdbool.h>

#include "nb/pde_bot/boundary_conditions/bcond.h"

typedef struct nb_bcond_iter_s nb_bcond_iter_t;

uint16_t nb_bcond_iter_get_memsize(void);
void nb_bcond_iter_init(void *iter_ptr);
void nb_bcond_iter_copy(void *iter_ptr, const void *src_iter_ptr);
void nb_bcond_iter_finish(void *iter_ptr);
void* nb_bcond_iter_create(void);
void* nb_bcond_iter_clone(const void *iter_ptr);
void nb_bcond_iter_destroy(void *iter_ptr);
void nb_bcond_iter_clear(void *iter_ptr);

void nb_bcond_iter_set_conditions(nb_bcond_iter_t *iter, nb_bcond_t *bcond,
				  nb_bcond_id type_id,
				  nb_bcond_where type_elem);
bool nb_bcond_iter_has_more(const nb_bcond_iter_t *iter);
void nb_bcond_iter_go_next(nb_bcond_iter_t *iter);

uint32_t nb_bcond_iter_get_id(const nb_bcond_iter_t *iter);
bool nb_bcond_iter_get_mask(const nb_bcond_iter_t *iter, uint8_t dof_id);
double nb_bcond_iter_get_val(const nb_bcond_iter_t *iter, uint8_t dof_id);

#endif
