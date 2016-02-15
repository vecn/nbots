#include <stdint.h>

#include "../container_struct.h"
#include "../iterator_struct.h"

#include "heap/avl_dst.h"
#include "heap/avl_iterator.h"

void heap_set_handlers(nb_container_t *container)
{
	container->b.init = heap_init;
	container->b.copy = heap_copy;
	container->b.clear = heap_clear;
	container->b.clone = heap_clone;
	container->b.destroy = heap_destroy;
	container->b.merge = heap_merge;
	container->b.insert = heap_insert;
	container->b.get_first = heap_get_first;
	container->b.delete_first = heap_delete_first;
	container->b.exist = heap_exist;
	container->b.delete = heap_delete;
	container->b.get_length = heap_get_length;
	container->b.is_empy = heap_is_empty;
	container->b.is_not_empty = heap_is_not_empty;
}

void heap_iterator_set_handlers(nb_iterator_t *iter)
{
	iter->b.init = heap_iter_init;
	iter->b.copy = heap_iter_copy;
	iter->b.clear = heap_iter_clear;
	iter->b.create = heap_iter_create;
	iter->b.clone = heap_iter_clone;
	iter->b.destroy = heap_iter_destroy;
	iter->b.restart = heap_iter_restart;
	iter->b.get_next = heap_iter_get_next;
	iter->b.has_more = heap_iter_has_more;
	iter->b.set_dst = heap_iter_set_dst;
}

void* heap_do(nb_container_t *container, const char* func,
		void *data, int8_t *status)
{
	*status = 2;
	return NULL;
}
