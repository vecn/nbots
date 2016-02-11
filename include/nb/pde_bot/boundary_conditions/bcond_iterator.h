#ifndef __NB_PDE_BOT_BCONDITIONS_ITERATOR_H__
#define __NB_PDE_BOT_BCONDITIONS_ITERATOR_H__

#include <stdint.h>
#include <stdbool.h>

#include "vcn/pde_bot/boundary_conditions.h"

typedef struct nb_bc_iterator_s nb_bc_iterator_t;

nb_bc_iterator_t* nb_bc_iterator_create(void);
void nb_bc_iterator_set_bcond(nb_bc_iterator_t *iter,
			      nb_bcond_t *bcond, nb_bcond_id type_id,
			      nb_bcond_where type_elem);
bool nb_bc_iterator_has_more(const nb_bc_iterator_t *iter);
uint32_t nb_bc_iterator_get_id(const nb_bc_iterator_t *iter);
bool nb_bc_iterator_view_mask(const nb_bc_iterator_t *iter, uint8_t dof_id);
double nb_bc_iterator_view_val(const nb_bc_iterator_t *iter, uint8_t dof_id);
void nb_bc_iterator_go_next(nb_bc_iterator_t *iter);
void nb_bc_iterator_destroy(nb_bc_iterator_t *iter);

#endif
