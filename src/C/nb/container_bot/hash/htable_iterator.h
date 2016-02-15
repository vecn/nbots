#ifndef __NB_HTABLE_ITERATOR_H__
#define __NB_HTABLE_ITERATOR_H__

#include <stdbool.h>

void* htable_iter_create(void);
void htable_iter_set_dst(void *iter_ptr, const void *const htable_ptr);
void* htable_iter_clone(const void *const iter_ptr);
void htable_iter_destroy(void *iter_ptr);
void htable_iter_restart(void *iter_ptr);

const void* htable_iter_get_next(void *iter_ptr);
bool htable_iter_has_more(const void *const iter_ptr);

#endif
