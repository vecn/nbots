#include <stdint.h>

#include "../container_struct.h"
#include "../iterator_struct.h"

#include "sorted/avl_dst.h"
#include "sorted/avl_iterator.h"

void sorted_set_handlers(nb_container_t *container)
{
	container->b.init = sorted_init;
	container->b.copy = sorted_copy;
	container->b.clear = sorted_clear;
	container->b.clone = sorted_clone;
	container->b.destroy = sorted_destroy;
	container->b.merge = sorted_merge;
	container->b.insert = sorted_insert;
	container->b.get_first = sorted_get_first;
	container->b.delete_first = sorted_delete_first;
	container->b.exist = sorted_exist;
	container->b.delete = sorted_delete;
	container->b.get_length = sorted_get_length;
	container->b.is_empy = sorted_is_empty;
	container->b.is_not_empty = sorted_is_not_empty;
}

void sorted_iterator_set_handlers(nb_iterator_t *iter)
{
	iter->b.init = sorted_iter_init;
	iter->b.copy = sorted_iter_copy;
	iter->b.clear = sorted_iter_clear;
	iter->b.create = sorted_iter_create;
	iter->b.clone = sorted_iter_clone;
	iter->b.destroy = sorted_iter_destroy;
	iter->b.restart = sorted_iter_restart;
	iter->b.get_next = sorted_iter_get_next;
	iter->b.has_more = sorted_iter_has_more;
	iter->b.set_dst = sorted_iter_set_dst;
}

void* sorted_do(nb_container_t *container, const char* func,
		void *data, int8_t *status)
{
	*status = 0;
	void *out = NULL;
	if (0 == strcmp("delete_last", func)) {
		out = avl_delete_last(container->cnt, container->op.key);
	} else {
		*status = 2;
	}
	return out;
}
