#include <stdint.h>

#include "../container_struct.h"
#include "../iterator_struct.h"

#include "hash/avl_dst.h"
#include "hash/avl_iterator.h"

void hash_set_handlers(nb_container_t *container)
{
	container->b.init = hash_init;
	container->b.copy = hash_copy;
	container->b.clear = hash_clear;
	container->b.clone = hash_clone;
	container->b.destroy = hash_destroy;
	container->b.merge = hash_merge;
	container->b.insert = hash_insert;
	container->b.get_first = hash_get_first;
	container->b.delete_first = hash_delete_first;
	container->b.exist = hash_exist;
	container->b.delete = hash_delete;
	container->b.get_length = hash_get_length;
	container->b.is_empy = hash_is_empty;
	container->b.is_not_empty = hash_is_not_empty;
}

void hash_iterator_set_handlers(nb_iterator_t *iter)
{
	iter->b.init = hash_iter_init;
	iter->b.copy = hash_iter_copy;
	iter->b.clear = hash_iter_clear;
	iter->b.create = hash_iter_create;
	iter->b.clone = hash_iter_clone;
	iter->b.destroy = hash_iter_destroy;
	iter->b.restart = hash_iter_restart;
	iter->b.get_next = hash_iter_get_next;
	iter->b.has_more = hash_iter_has_more;
	iter->b.set_dst = hash_iter_set_dst;
}

void* hash_do(nb_container_t *container, const char* func,
		void *data, int8_t *status)
{
	*status = 0;
	void *out = NULL;
	if (0 == strcmp("get_size", func)) {
		uint32_t *size = malloc(sizeof(*size));
		*size = htable_get_size(container->cnt);
		out = size;
	} else if (0 == strcmp("get_N_collisions", func)) {
		uint32_t *N = malloc(sizeof(*N));
		*N = htable_get_N_collisions(container->cnt);
		out = N;
	} else if (0 == strcmp("get_collisions", func)) {
		out = htable_get_collisions(container->cnt);
	} else {
		*status = 2;
	}
	return out;
}
