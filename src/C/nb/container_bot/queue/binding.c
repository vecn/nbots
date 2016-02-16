#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "../container_struct.h"
#include "../iterator_struct.h"

#include "queue_dst.h"
#include "queue_iterator.h"
#include "binding.h"

void queue_set_handlers(nb_container_t *container)
{
	container->b.init = queue_init;
	container->b.copy = queue_copy;
	container->b.clear = queue_clear;
	container->b.clone = queue_clone;
	container->b.destroy = queue_destroy;
	container->b.merge = queue_merge;
	container->b.insert = queue_insert;
	container->b.get_first = queue_get_first;
	container->b.delete_first = queue_delete_first;
	container->b.exist = queue_exist;
	container->b.delete = queue_delete;
	container->b.get_length = queue_get_length;
	container->b.is_empty = queue_is_empty;
	container->b.is_not_empty = queue_is_not_empty;
}

void queue_iterator_set_handlers(nb_iterator_t *iter)
{
	iter->b.init = queue_iter_init;
	iter->b.copy = queue_iter_copy;
	iter->b.clear = queue_iter_clear;
	iter->b.create = queue_iter_create;
	iter->b.clone = queue_iter_clone;
	iter->b.destroy = queue_iter_destroy;
	iter->b.restart = queue_iter_restart;
	iter->b.get_next = queue_iter_get_next;
	iter->b.has_more = queue_iter_has_more;
	iter->b.set_dst = queue_iter_set_dst;
}

void* queue_do(nb_container_t *container, const char* func,
	       void *data, int8_t *status)
{
	*status = 0;
	void *out = NULL;
	if (0 == strcmp("insert_first", func)) {
		queue_insert_first(container->cnt, data, container->op.key);
	} else {
		*status = 2;
	}
	return out;
}
