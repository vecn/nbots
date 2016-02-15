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

inline void* stack_iter_create(void)
{
	return calloc(1, sizeof(iter_t));
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

void* stack_iter_clone(const void *const iter_ptr)
{
	const iter_t *const restrict iter = iter_ptr;
	iter_t *const restrict clone = calloc(1, sizeof(*clone));
	clone->start = iter->start;
	clone->node = iter->node;
	clone->is_init = iter->is_init;
	return clone;
}

inline void stack_iter_destroy(void *iter_ptr)
{
	free(iter_ptr);
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
