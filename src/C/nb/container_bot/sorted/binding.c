#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "../container_struct.h"
#include "../iterator_struct.h"

#include "avl_dst.h"
#include "avl_iterator.h"
#include "binding.h"

void sorted_set_handlers(nb_container_t *container)
{
	container->b.init = avl_init;
	container->b.copy = avl_copy;
	container->b.finish = avl_finish;
	container->b.clone = avl_clone;
	container->b.destroy = avl_destroy;
	container->b.clear = avl_clear;
	container->b.merge = avl_merge;
	container->b.insert = avl_insert;
	container->b.get_first = avl_get_first;
	container->b.delete_first = avl_delete_first;
	container->b.exist = avl_exist;
	container->b.delete = avl_delete;
	container->b.get_length = avl_get_length;
	container->b.is_empty = avl_is_empty;
	container->b.is_not_empty = avl_is_not_empty;
}

void sorted_iterator_set_handlers(nb_iterator_t *iter)
{
	iter->b.init = avl_iter_init;
	iter->b.copy = avl_iter_copy;
	iter->b.finish = avl_iter_finish;
	iter->b.create = avl_iter_create;
	iter->b.clone = avl_iter_clone;
	iter->b.destroy = avl_iter_destroy;
	iter->b.clear = avl_iter_clear;
	iter->b.restart = avl_iter_restart;
	iter->b.get_next = avl_iter_get_next;
	iter->b.has_more = avl_iter_has_more;
	iter->b.set_dst = avl_iter_set_dst;
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
