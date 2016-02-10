#ifndef __VCN_AVL_ITERATOR_H__
#define __VCN_AVL_ITERATOR_H__

#include <stdbool.h>

void* avl_iter_create(void);
void avl_iter_set_dst(void *iter_ptr, const void *const avl_ptr);
void* avl_iter_clone(const void *const iter_ptr);
void avl_iter_destroy(void *iter_ptr);
void avl_iter_restart(void *iter_ptr);

const void* avl_iter_get_next(void *iter_ptr);
bool avl_iter_has_more(const void *const iter_ptr);

#endif
