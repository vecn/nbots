#include <stdint.h>

#include "../container_struct.h"
#include "../iterator_struct.h"

#include "stack/stack_dst.h"
#include "stack/stack_iterator.h"

void stack_set_handlers(nb_container_t *container)
{
	container->b.init = stack_init;
	container->b.copy = stack_copy;
	container->b.clear = stack_clear;
	container->b.clone = stack_clone;
	container->b.destroy = stack_destroy;
	container->b.merge = stack_merge;
	container->b.insert = stack_insert;
	container->b.get_first = stack_get_first;
	container->b.delete_first = stack_delete_first;
	container->b.exist = stack_exist;
	container->b.delete = stack_delete;
	container->b.get_length = stack_get_length;
	container->b.is_empy = stack_is_empty;
	container->b.is_not_empty = stack_is_not_empty;
}

void stack_iterator_set_handlers(nb_iterator_t *iter)
{
	iter->b.init = stack_iter_init;
	iter->b.copy = stack_iter_copy;
	iter->b.clear = stack_iter_clear;
	iter->b.create = stack_iter_create;
	iter->b.clone = stack_iter_clone;
	iter->b.destroy = stack_iter_destroy;
	iter->b.restart = stack_iter_restart;
	iter->b.get_next = stack_iter_get_next;
	iter->b.has_more = stack_iter_has_more;
	iter->b.set_dst = stack_iter_set_dst;
}

void* stack_do(nb_container_t *container, const char* func,
		void *data, int8_t *status)
{
	*status = 2
	return NULL;
}
