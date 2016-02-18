/******************************************************************************
 *   List DST: Circular Linked List                                           *
 *   2011-2015 Victor Eduardo Cardoso Nungaray                                *
 *   Twitter: @victore_cardoso                                                *
 *   email: victorc@cimat.mx                                                  *
 ******************************************************************************/

#include <stdlib.h>
#include <string.h>

#include "stack_node.h"
#include "stack_struct.h"
#include "stack_dst.h"

static bool is_not_empty(const stack_t *const list);
static node_t* get_first(const stack_t *const list);
static void* malloc_stack(void);
static void insert_node_as_starting(stack_t *list, const void *const val);
static void add_node(stack_t *list, node_t *node);
static void add_first_node(stack_t *list, node_t *node);
static void null_destroy(void *val);
static void link_node(stack_t *stack, node_t *node);
static node_t* exist_node(const stack_t *const stack, const void *val,
			  int8_t (*compare)(const void*, const void*));
static void unlink_node(stack_t *stack, const node_t *const node);

inline uint16_t stack_get_memsize(void)
{
	return sizeof(stack_t);
}

inline void stack_init(void *stack_ptr)
{
	stack_t *stack = stack_ptr;
	stack->length = 0;
	stack->end = NULL;
}

void stack_copy(void *stack_ptr, const void *src_stack_ptr,
		void* (*clone)(const void*))
{
	stack_t *stack = stack_ptr;
	const stack_t *src_stack = src_stack_ptr;
	
	stack->length = src_stack->length;

	if (is_not_empty(src_stack)) {
		node_t *end = node_clone(src_stack->end, clone);
		stack->end = end;

		node_t *iter = get_first(src_stack);
		do {
			node_t *node = node_clone(iter, clone);
			end->next = node;
			end = node;
			iter = iter->next;
		} while (iter != src_stack->end);
		end->next = stack->end;
	} else {
		stack->end = NULL;
	}
}


static inline bool is_not_empty(const stack_t *const restrict stack)
{
	return (stack->end != NULL);
}

static inline node_t* get_first(const stack_t *const restrict stack)
{
	return stack->end->next;
}

inline void stack_finish(void *stack_ptr,
			 void (*destroy)(void*))
{
	stack_clear(stack_ptr, destroy);
}

inline void* stack_create(void)
{
	void *stack = malloc_stack();
	stack_init(stack);
	return stack;
}

static inline void* malloc_stack(void)
{
	uint16_t size = stack_get_memsize();
  	return malloc(size);
}

inline void* stack_clone(const void *const stack_ptr,
			 void* (*clone)(const void*))
{
	void *stack = malloc_stack();
	stack_copy(stack, stack_ptr, clone);
	return stack;
}

inline void stack_destroy(void *stack_ptr,
		  void (*destroy)(void*))
{
	stack_finish(stack_ptr, destroy);
	free(stack_ptr);
}

void stack_clear(void *stack_ptr,
		 void (*destroy)(void*))
{
	stack_t *stack = stack_ptr;
	if (is_not_empty(stack)) {
		node_t* iter = get_first(stack);
		stack->end->next = NULL;
		while (NULL != iter) {
			node_t* to_destroy = iter;
			iter = iter->next;
			node_destroy(to_destroy, destroy);
		}
		stack->length = 0;
		stack->end = NULL;
	}
}

void stack_merge(void *stack1_ptr, void *stack2_ptr,
		 uint32_t (*key)(const void*),
		 int8_t (*compare)(const void*, const void*))
{
	stack_t *stack1 = stack1_ptr;
	stack_t *stack2 = stack2_ptr;
	if (is_not_empty(stack2)) {
		stack1->length += stack2->length;
		if (is_not_empty(stack1)) {
			node_t *first_node = get_first(stack1);
			stack1->end->next = get_first(stack2);
			stack2->end->next = first_node;
		}
		stack1->end = stack2->end;
		stack2->length = 0;
		stack2->end = NULL;
	}
}

inline bool stack_insert(void *stack_ptr, const void *val,
			 uint32_t (*key)(const void*),
			 int8_t (*compare)(const void*, const void*))
{
	stack_t *stack = stack_ptr;
	insert_node_as_starting(stack, val);
	return true;
}

static inline void insert_node_as_starting(stack_t *restrict stack,
					   const void *const restrict val)
{
	node_t *const restrict node = node_create();
	node->val = (void*) val;
	add_node(stack, node);
	stack->length += 1;
}

static inline void add_node(stack_t *restrict stack,
			    node_t *restrict node)
{
	if (stack_is_empty(stack))
		add_first_node(stack, node);
	else
		link_node(stack, node);
}

static inline void add_first_node(stack_t *restrict stack,
				  node_t *restrict node)
{
	stack->end = node;
	node->next = node;
}

static inline void link_node(stack_t *restrict stack, 
			     node_t *restrict node)
{
	node->next = get_first(stack);
	stack->end->next = node;
}

inline void* stack_get_first(const void *const stack_ptr)
{
	const stack_t *const stack = stack_ptr;
	void *val = NULL;
	if (is_not_empty(stack)) {
		node_t *first = get_first(stack);
		val = first->val;
	}
	return val;
}

void* stack_delete_first(void *stack_ptr,
			uint32_t (*key)(const void*))
{
	stack_t *stack = stack_ptr;
	void *val = NULL;
	if (is_not_empty(stack)) {
		node_t *first = get_first(stack);
		val = first->val;
		if (stack->end != first)
		  stack->end->next = first->next;
		else
		  stack->end = NULL;
		node_destroy(first, null_destroy);
		stack->length -= 1;
	}
	return val;
}

static inline void null_destroy(void *val)
{
  ; /* Null statement */
}

void* stack_exist(const void *const stack_ptr, const void *val,
		  uint32_t (*key)(const void*),
		  int8_t (*compare)(const void*, const void*))
{
	const stack_t *const stack = stack_ptr;
	void *existing_val = NULL;
	if (is_not_empty(stack)) {
		node_t *node = exist_node(stack, val, compare);
		if (NULL != node)
			existing_val = node->val;
	}
	return existing_val;
}

static node_t* exist_node(const stack_t *const stack, const void *val,
			  int8_t (*compare)(const void*, const void*))
{
	node_t *first = get_first(stack);
	node_t *existing_node = NULL;
	if (0 == compare(first->val, val)) {
		existing_node = first;
	} else {
		node_t *iter = first->next;
		while (iter != first) {
			if (0 == compare(iter->val, val)) {
				existing_node = iter;
				break;
			}
			iter = iter->next;
		}
	}
	return existing_node;
}

void* stack_delete(void *stack_ptr, const void *val,
		  uint32_t (*key)(const void*),
		  int8_t (*compare)(const void*, const void*))
{
	stack_t *stack = stack_ptr;
	void *deleted_val = NULL;
	if (is_not_empty(stack)) {
		node_t *node = exist_node(stack, val, compare);
		if (NULL != node) {
			unlink_node(stack, node);
			deleted_val = node->val;
			node_destroy(node, null_destroy);
			stack->length -= 1;
		}
	}
	return deleted_val;
}

static void unlink_node(stack_t *stack, const node_t *const node)
{
	node_t *prev = node_get_prev(node);
	if (node == prev /* implies node is equal to stack->end */)
		stack->end = NULL;
	else if (node == stack->end /* implies node is different to prev */)
		stack->end = prev;
	prev->next = node->next;
}

inline uint32_t stack_get_length(const void *const stack_ptr)
{
	const stack_t *const stack = stack_ptr;
	return stack->length;
}

inline bool stack_is_empty(const void *const stack_ptr)
{
	const stack_t *const stack = stack_ptr;
	return (NULL == stack->end);
}

inline bool stack_is_not_empty(const void *const stack_ptr)
{
	return is_not_empty(stack_ptr);
}
