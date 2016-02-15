#ifndef __NB_CONTAINER_BOT_STACK_STACK_ITERATOR_H__
#define __NB_CONTAINER_BOT_STACK_STACK_ITERATOR_H__

#include <stdbool.h>
#include <stdint.h>

uint16_t stack_get_iter_memsize(void);
void stack_iter_init(void *iter_ptr);
void stack_iter_copy(void *iter_ptr, const void *srciter_ptr);
void stack_iter_clear(void *iter_ptr);
void* stack_iter_create(void);
void* stack_iter_clone(const void *const iter_ptr);
void stack_iter_destroy(void *iter_ptr);

void stack_iter_set_dst(void *iter_ptr, const void *const stack_ptr);
void stack_iter_restart(void *iter_ptr);

const void* stack_iter_get_next(void *iter_ptr);
bool stack_iter_has_more(const void *const iter_ptr);

#endif
