#ifndef __NB_CONTAINER_BOT_QUEUE_NB_QUEUE_ITERATOR_H__
#define __NB_CONTAINER_BOT_QUEUE_NB_QUEUE_ITERATOR_H__

#include <stdbool.h>
#include <stdint.h>

uint16_t nb_queue_iter_get_memsize(void);
void nb_queue_iter_init(void *iter_ptr);
void nb_queue_iter_copy(void *iter_ptr, const void *src_iter_ptr);
void nb_queue_iter_finish(void *iter_ptr);
void* nb_queue_iter_create(void);
void* nb_queue_iter_clone(const void *const iter_ptr);
void nb_queue_iter_destroy(void *iter_ptr);
void nb_queue_iter_clear(void *iter_ptr);

void nb_queue_iter_set_dst(void *iter_ptr, const void *const queue_ptr);
void nb_queue_iter_restart(void *iter_ptr);

const void* nb_queue_iter_get_next(void *iter_ptr);
bool nb_queue_iter_has_more(const void *const iter_ptr);

#endif
