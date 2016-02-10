#ifndef __NB_LIST_ITERATOR_H__
#define __NB_LIST_ITERATOR_H__

#include <stdbool.h>

void* list_iter_create(void);
void list_iter_set_dst(void *iter_ptr, const void *const list_ptr);
void* list_iter_clone(const void *const iter_ptr);
void list_iter_destroy(void *iter_ptr);
void list_iter_restart(void *iter_ptr);

const void* list_iter_get_next(void *iter_ptr);
bool list_iter_has_more(const void *const iter_ptr);

#endif
