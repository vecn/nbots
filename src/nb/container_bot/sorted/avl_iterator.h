#ifndef __NB_CONTAINER_BOT_SORTED_AVL_ITERATOR_H__
#define __NB_CONTAINER_BOT_SORTED_AVL_ITERATOR_H__

#include <stdint.h>
#include <stdbool.h>

uint16_t avl_iter_get_memsize(void);
void avl_iter_init(void *iter_ptr);
void avl_iter_copy(void *iter_ptr, const void *src_iter_ptr);
void avl_iter_finish(void *iter_ptr);
void* avl_iter_create(void);
void* avl_iter_clone(const void *iter_ptr);
void avl_iter_destroy(void *iter_ptr);
void avl_iter_clear(void *iter_ptr);

void avl_iter_set_dst(void *iter_ptr, const void *avl_ptr);
void avl_iter_restart(void *iter_ptr);

const void* avl_iter_get_next(void *iter_ptr);
bool avl_iter_has_more(const void *const iter_ptr);

#endif
