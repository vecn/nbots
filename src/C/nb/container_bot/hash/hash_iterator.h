#ifndef __CONTAINER_BOT_HASH_HASH_ITERATOR_H__
#define __CONTAINER_BOT_HASH_HASH_ITERATOR_H__

#include <stdbool.h>

uint16_t hash_iter_get_memsize(void);
void hash_iter_init(void *iter_ptr);
void hash_iter_copy(void *iter_ptr, const void *src_iter_ptr);
void hash_iter_finish(void *iter_ptr);
void* hash_iter_create(void);
void* hash_iter_clone(const void *iter_ptr);
void hash_iter_destroy(void *iter_ptr);
void hash_iter_clear(void *iter_ptr);

void hash_iter_set_dst(void *iter_ptr, const void *hash_ptr);
void hash_iter_restart(void *iter_ptr);

const void* hash_iter_get_next(void *iter_ptr);
bool hash_iter_has_more(const void *const iter_ptr);

#endif
