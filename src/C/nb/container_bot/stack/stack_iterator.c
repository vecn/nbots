/******************************************************************************
 *   Stack iterator                                                           *
 *   2011-2015 Victor Eduardo Cardoso Nungaray                                *
 *   Twitter: @victore_cardoso                                                *
 *   email: victorc@cimat.mx                                                  *
 ******************************************************************************/

#include <stdlib.h>
#include <string.h>

#include "stack_node.h"
#include "stack_struct.h"
#include "stack_iterator.h"

typedef struct {
	bool is_init;
	node_t *node;
	const node_t *start;
} iter_t;

static void* malloc_iter(void);

inline uint16_t stack_iter_get_memsize(void)
{
	return sizeof(iter_t);
}

inline void stack_iter_init(void *iter_ptr)
{
	memset(iter_ptr, 0, stack_iter_get_memsize());
}

void stack_iter_copy(void *iter_ptr, const void *src_iter_ptr)
{
	iter_t *iter = iter_ptr;
	const iter_t *src_iter = src_iter_ptr;
	iter->start = src_iter->start;
	iter->node = src_iter->node;
	iter->is_init = src_iter->is_init;
}

inline void stack_iter_clear(void *iter_ptr)
{
	memset(iter_ptr, 0, stack_iter_get_memsize());
}

inline void* stack_iter_create(void)
{
	void *iter = malloc_iter();
	stack_iter_init(iter);
	return iter;
}

static inline void* malloc_iter(void)
{
	uint16_t size = stack_iter_get_memsize();
	return malloc(size);
}

inline void* stack_iter_clone(const void *const iter_ptr)
{
	void *iter = malloc_iter();
	stack_iter_copy(iter, iter_ptr);
	return iter;
}

inline void stack_iter_destroy(void *iter_ptr)
{
	free(iter_ptr);
}

void stack_iter_set_dst(void *iter_ptr, const void *const stack_ptr)
{
	iter_t *iter = iter_ptr;
	if (NULL != stack_ptr) {
		const stack_t* stack = stack_ptr;
		iter->start = stack->end->next;
	} else {
		iter->start = NULL;
	}
	stack_iter_restart(iter);
}

inline void stack_iter_restart(void* iter_ptr)
{
	iter_t *iter = iter_ptr;
	iter->node = (node_t*) iter->start;
	iter->is_init = (NULL != iter->start);
}

const void* stack_iter_get_next(void *iter_ptr)
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

inline bool stack_iter_has_more(const void *const iter_ptr)
{
	const iter_t *const restrict iter = iter_ptr;
	return (iter->start != iter->node) || iter->is_init;
}
