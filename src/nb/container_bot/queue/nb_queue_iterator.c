#include <stdlib.h>
#include <string.h>

#include "nb/memory_bot.h"

#include "nb_queue_node.h"
#include "nb_queue_struct.h"
#include "nb_queue_iterator.h"

typedef struct {
	bool is_init;
	nb_queue_node_t *node;
	const nb_queue_node_t *start;
} iter_t;

static void* allocate_iter(void);

uint16_t nb_queue_iter_get_memsize(void)
{
	return sizeof(iter_t);
}

void nb_queue_iter_init(void *iter_ptr)
{
	memset(iter_ptr, 0, nb_queue_iter_get_memsize());
}

void nb_queue_iter_copy(void *iter_ptr, const void *src_iter_ptr)
{
	iter_t *iter = iter_ptr;
	const iter_t *src_iter = src_iter_ptr;
	iter->start = src_iter->start;
	iter->node = src_iter->node;
	iter->is_init = src_iter->is_init;
}

void nb_queue_iter_finish(void *iter_ptr)
{
	nb_queue_iter_clear(iter_ptr);
}

void* nb_queue_iter_create(void)
{
	void *iter = allocate_iter();
	nb_queue_iter_init(iter);
	return iter;
}

static inline void* allocate_iter(void)
{
	uint16_t size = nb_queue_iter_get_memsize();
	return nb_allocate_mem(size);
}

void* nb_queue_iter_clone(const void *const iter_ptr)
{
	void *iter = allocate_iter();
	nb_queue_iter_copy(iter, iter_ptr);
	return iter;
}

void nb_queue_iter_destroy(void *iter_ptr)
{
	nb_queue_iter_finish(iter_ptr);
	nb_free_mem(iter_ptr);
}

void nb_queue_iter_clear(void *iter_ptr)
{
	memset(iter_ptr, 0, nb_queue_iter_get_memsize());
}

void nb_queue_iter_set_dst(void *iter_ptr, const void *const queue_ptr)
{
	iter_t *iter = iter_ptr;
	iter->start = NULL;
	if (NULL != queue_ptr) {
		const nb_queue_t* queue = queue_ptr;
		if (NULL != queue->end)
			iter->start = queue->end->next;
	}
	nb_queue_iter_restart(iter);
}


void nb_queue_iter_restart(void* iter_ptr)
{
	iter_t *iter = iter_ptr;
	iter->node = (nb_queue_node_t*) iter->start;
	iter->is_init = (NULL != iter->start);
}

const void* nb_queue_iter_get_next(void *iter_ptr)
{
	iter_t *iter = iter_ptr;
	void *val = NULL;
	if (NULL != iter->node) {
		val = iter->node->val;
		iter->node = iter->node->next;
		iter->is_init = false;
	}
	return val;
}

bool nb_queue_iter_has_more(const void *const iter_ptr)
{
	const iter_t *iter = iter_ptr;
	return (iter->start != iter->node) || iter->is_init;
}
