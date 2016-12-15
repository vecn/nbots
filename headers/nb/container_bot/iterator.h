/**
 * @file iterator.h
 * @brief  Generic iterator.
 * @author Victor Eduardo Cardoso Nungaray
 * @n victorc@@cimat.mx
 * @n @@victore_cardoso
 */

#ifndef __NB_CONTAINER_BOT_ITERATOR_H__
#define __NB_CONTAINER_BOT_ITERATOR_H__

#include <stdbool.h>
#include <stdint.h>
#include "nb/container_bot/container.h"
  
typedef struct nb_iterator_s nb_iterator_t;

uint16_t nb_iterator_get_memsize(void);
void nb_iterator_init(void *iter_ptr);
void nb_iterator_copy(void *iter_ptr, const void *src_iter_ptr);
void nb_iterator_finish(void* iter_ptr);
void* nb_iterator_create(void);
void* nb_iterator_clone(const void *iter_ptr);
void nb_iterator_destroy(void *iter_ptr);
void nb_iterator_clear(void* iter_ptr);

void nb_iterator_set_container(nb_iterator_t *iter,
			       const nb_container_t *const container);
void nb_iterator_restart(nb_iterator_t *iter);
const void* nb_iterator_get_next(nb_iterator_t *iter);
bool nb_iterator_has_more(const nb_iterator_t *const iter);

#endif
