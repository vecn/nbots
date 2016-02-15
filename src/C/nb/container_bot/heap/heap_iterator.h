#ifndef __NB_NB_CONTAINER_HEAP_HEAP_ITERATOR_H__
#define __NB_NB_CONTAINER_HEAP_HEAP_ITERATOR_H__

#include <stdint.h>
#include <stdbool.h>

uint16_t heap_iter_get_memsize(void);
void heap_iter_init(void *iter_ptr);
void heap_iter_copy(void *iter_ptr, const void *src_iter_ptr);
void heap_iter_clear(void *iter_ptr);
void* heap_iter_create(void);
void* heap_iter_clone(const void *const iter_ptr);
void heap_iter_destroy(void *iter_ptr);

void heap_iter_set_dst(void *iter_ptr, const void *const heap_ptr);
void heap_iter_restart(void *iter_ptr);

const void* heap_iter_get_next(void *iter_ptr);
bool heap_iter_has_more(const void *const iter_ptr);

#endif
