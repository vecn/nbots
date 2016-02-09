/******************************************************************************
 *   Heap iterator                                                            *
 *   2011-2015 Victor Eduardo Cardoso Nungaray                                *
 *   Twitter: @victore_cardoso                                                *
 *   email: victorc@cimat.mx                                                  *
 ******************************************************************************/

#include <stdlib.h>
#include <string.h>

#include "heap_tree.h"
#include "heap_dst.h"
#include "heap_iterator.h"

typedef struct {
	/* Pre-order iterator */
	const htree_t *start;
	htree_t *prev;
	htree_t *current;
} iter_t;

static void go_up_while_down_is_visited(iter_t *iter);
static bool its_possible_going_down(const iter_t *const iter);
static void go_down_if_not_visited(iter_t *iter);
static bool down_is_visited(const iter_t *const iter);
static void go_down(iter_t *iter);
static void go_right(iter_t *iter);
static void go_up_may_be_right(iter_t *iter);

inline void* heap_iter_create(void)
{
	return calloc(1, sizeof(iter_t));
}
 
void heap_iter_set_dst(void *iter_ptr, const void *const  heap_ptr)
{
	iter_t *iter = iter_ptr;
	if (NULL != heap_ptr)
		iter->start = heap_get_iterator_start(heap_ptr);
	heap_iter_restart(iter);
}

void* heap_iter_clone(const void *const iter_ptr)
{
	const iter_t *const iter = iter_ptr;
	iter_t *clone = heap_iter_create();
	clone->start = iter->start;
	clone->prev = iter->prev;
	clone->current = iter->current;
	return clone;
}

inline void heap_iter_destroy(void *iter_ptr)
{
	free(iter_ptr);
}

inline void heap_iter_restart(void *iter_ptr)
{
	iter_t *iter = iter_ptr;
	iter->prev = NULL;
	iter->current = (htree_t*) iter->start;
}

const void* heap_iter_get_next(void *iter_ptr)
{
	iter_t *iter = iter_ptr;
	void *val = NULL;
	if (NULL != iter->current) {
		val = iter->current->val;
		if (its_possible_going_down(iter)) {
			go_down_if_not_visited(iter);
		} else {
			go_right(iter);
			go_up_while_down_is_visited(iter);
		}
	}
	return val;
}

static inline void go_up_while_down_is_visited(iter_t *iter)
{
	while (down_is_visited(iter))
		go_up_may_be_right(iter);
}

static inline bool its_possible_going_down(const iter_t *const iter)
{
	return (NULL != iter->current->left);
}

static void go_down_if_not_visited(iter_t *iter)
{
	if (down_is_visited(iter))
		go_up_while_down_is_visited(iter);
	else
		go_down(iter);
}

static bool down_is_visited(const iter_t *const iter)
{
	bool last_was_up = false;
	bool prev_and_current_not_null = 
		(NULL != iter->prev) && (NULL != iter->current);
	if (prev_and_current_not_null) {
		bool left_not_null = (NULL != iter->current->left);
		if (left_not_null) {
			htree_t* left = iter->current->left;
			bool prev_is_ltree = (left == iter->prev);
			bool prev_is_rtree = (left->right == iter->prev);
			last_was_up = prev_is_ltree || prev_is_rtree;
		}
	}
	return last_was_up;
}

static inline void go_down(iter_t *restrict iter)
{
	iter->prev = iter->current;
	iter->current = iter->current->left;
}

static inline void go_right(iter_t *restrict iter)
{
	iter->prev = iter->current;
	iter->current = iter->current->right;
}

static inline void go_up_may_be_right(iter_t *restrict iter)
{
	go_right(iter);
}

inline bool heap_iter_has_more(const void *const iter_ptr)
{
	const iter_t *const restrict iter = iter_ptr;
	return (NULL != iter->current);
}
