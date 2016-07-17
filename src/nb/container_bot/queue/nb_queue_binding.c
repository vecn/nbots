#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "../container_struct.h"
#include "../iterator_struct.h"

#include "nb_queue_dst.h"
#include "nb_queue_iterator.h"
#include "nb_queue_binding.h"

void nb_queue_set_handlers(nb_container_t *container)
{
	container->b.init = nb_queue_init;
	container->b.copy = nb_queue_copy;
	container->b.finish = nb_queue_finish;
	container->b.clone = nb_queue_clone;
	container->b.destroy = nb_queue_destroy;
	container->b.clear = nb_queue_clear;
	container->b.merge = nb_queue_merge;
	container->b.insert = nb_queue_insert;
	container->b.get_first = nb_queue_get_first;
	container->b.delete_first = nb_queue_delete_first;
	container->b.exist = nb_queue_exist;
	container->b.delete = nb_queue_delete;
	container->b.get_length = nb_queue_get_length;
	container->b.is_empty = nb_queue_is_empty;
	container->b.is_not_empty = nb_queue_is_not_empty;
}

void nb_queue_iterator_set_handlers(nb_iterator_t *iter)
{
	iter->b.init = nb_queue_iter_init;
	iter->b.copy = nb_queue_iter_copy;
	iter->b.finish = nb_queue_iter_finish;
	iter->b.create = nb_queue_iter_create;
	iter->b.clone = nb_queue_iter_clone;
	iter->b.destroy = nb_queue_iter_destroy;
	iter->b.clear = nb_queue_iter_clear;
	iter->b.restart = nb_queue_iter_restart;
	iter->b.get_next = nb_queue_iter_get_next;
	iter->b.has_more = nb_queue_iter_has_more;
	iter->b.set_dst = nb_queue_iter_set_dst;
}

void* nb_queue_do(nb_container_t *container, const char* func,
		  void *data, int8_t *status)
{
	*status = 0;
	void *out = NULL;
	if (0 == strcmp("insert_first", func)) {
		nb_queue_insert_first(container->cnt, data,
				      container->op.key,
				      container->op.compare);
	} else {
		*status = 2;
	}
	return out;
}
