/******************************************************************************
 *   List DST: Circular Linked List                                           *
 *   2011-2015 Victor Eduardo Cardoso Nungaray                                *
 *   Twitter: @victore_cardoso                                                *
 *   email: victorc@cimat.mx                                                  *
 ******************************************************************************/

#include <stdlib.h>
#include <string.h>

#include "nb/memory_bot.h"

#include "stack_node.h"
#include "stack_struct.h"
#include "stack_dst.h"

static bool is_not_empty(const nb_stack_t *const list);
static stack_node_t* get_first(const nb_stack_t *const list);
static void* malloc_stack(void);
static void destroy_values(nb_stack_t *stack, void (*destroy)(void*));
static void insert_stack_node_as_starting(nb_stack_t *list,
					  const void *const val);
static void add_node(nb_stack_t *list, stack_node_t *node);
static void add_first_node(nb_stack_t *list, stack_node_t *node);
static void link_node(nb_stack_t *stack, stack_node_t *node);
static stack_node_t* exist_node(const nb_stack_t *const stack, const void *val,
			  int8_t (*compare)(const void*, const void*));
static void unlink_node(nb_stack_t *stack, const stack_node_t *const node);

uint16_t stack_get_memsize(void)
{
	return sizeof(nb_stack_t) + nb_membank_get_memsize();
}

void stack_init(void *stack_ptr)
{
	char *memblock = stack_ptr;
	nb_stack_t *stack = (void*) memblock;
	memset(stack_ptr, 0, stack_get_memsize());
	stack->length = 0;
	stack->end = NULL;
	stack->membank = (void*) (memblock  + sizeof(nb_stack_t));
	nb_membank_init(stack->membank, stack_node_get_memsize());
}

void stack_copy(void *stack_ptr, const void *src_stack_ptr,
		void* (*clone)(const void*))
{
	nb_stack_t *stack = stack_ptr;
	const nb_stack_t *src_stack = src_stack_ptr;

	stack_init(stack);
	
	stack->length = src_stack->length;

	if (is_not_empty(src_stack)) {
		stack_node_t *end = stack_node_clone(stack->membank,
						     src_stack->end, clone);
		stack->end = end;

		stack_node_t *iter = get_first(src_stack);
		do {
			stack_node_t *node = stack_node_clone(stack->membank,
							      iter, clone);
			end->next = node;
			end = node;
			iter = iter->next;
		} while (iter != src_stack->end);
		end->next = stack->end;
	} else {
		stack->end = NULL;
	}
}


static inline bool is_not_empty(const nb_stack_t *const restrict stack)
{
	return (stack->end != NULL);
}

static inline stack_node_t* get_first(const nb_stack_t *const restrict stack)
{
	return stack->end->next;
}

void stack_finish(void *stack_ptr, void (*destroy)(void*))
{
	nb_stack_t *stack = stack_ptr;
	destroy_values(stack, destroy);
	nb_membank_finish(stack->membank);
}

void* stack_create(void)
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

void* stack_clone(const void *const stack_ptr,
		  void* (*clone)(const void*))
{
	void *stack = malloc_stack();
	stack_copy(stack, stack_ptr, clone);
	return stack;
}

void stack_destroy(void *stack_ptr, void (*destroy)(void*))
{
	stack_finish(stack_ptr, destroy);
	free(stack_ptr);
}

void stack_clear(void *stack_ptr, void (*destroy)(void*))
{
	nb_stack_t *stack = stack_ptr;
	destroy_values(stack, destroy);
	nb_membank_clear(stack->membank);
	stack->length = 0;
	stack->end = NULL;
}

static void destroy_values(nb_stack_t *stack, void (*destroy)(void*))
{
	if (NULL != destroy) {
		if (is_not_empty(stack)) {
			stack_node_t* iter = get_first(stack);
			stack->end->next = NULL;
			while (NULL != iter) {
				destroy(iter->val);
				iter = iter->next;
			}
		}
	}
}

void stack_merge(void *stack1_ptr, void *stack2_ptr,
		 uint32_t (*key)(const void*),
		 int8_t (*compare)(const void*, const void*))
{
	nb_stack_t *stack1 = stack1_ptr;
	nb_stack_t *stack2 = stack2_ptr;
	if (is_not_empty(stack2)) {
		stack1->length += stack2->length;
		if (is_not_empty(stack1)) {
			stack_node_t *first_node = get_first(stack1);
			stack1->end->next = get_first(stack2);
			stack2->end->next = first_node;
		}
		stack1->end = stack2->end;
		stack2->length = 0;
		stack2->end = NULL;

		nb_membank_merge(stack1->membank, stack2->membank);
	}
}

bool stack_insert(void *stack_ptr, const void *val,
		  uint32_t (*key)(const void*),
		  int8_t (*compare)(const void*, const void*))
{
	nb_stack_t *stack = stack_ptr;
	insert_stack_node_as_starting(stack, val);
	return true;
}

static inline void insert_stack_node_as_starting(nb_stack_t * stack,
						 const void *const val)
{
	stack_node_t *const restrict node = stack_node_create(stack->membank);
	node->val = (void*) val;
	add_node(stack, node);
	stack->length += 1;
}

static inline void add_node(nb_stack_t *restrict stack,
			    stack_node_t *restrict node)
{
	if (stack_is_empty(stack))
		add_first_node(stack, node);
	else
		link_node(stack, node);
}

static inline void add_first_node(nb_stack_t *restrict stack,
				  stack_node_t *restrict node)
{
	stack->end = node;
	node->next = node;
}

static inline void link_node(nb_stack_t *restrict stack, 
			     stack_node_t *restrict node)
{
	node->next = get_first(stack);
	stack->end->next = node;
}

void* stack_get_first(const void *const stack_ptr)
{
	const nb_stack_t *const stack = stack_ptr;
	void *val = NULL;
	if (is_not_empty(stack)) {
		stack_node_t *first = get_first(stack);
		val = first->val;
	}
	return val;
}

void* stack_delete_first(void *stack_ptr,
			uint32_t (*key)(const void*))
{
	nb_stack_t *stack = stack_ptr;
	void *val = NULL;
	if (is_not_empty(stack)) {
		stack_node_t *first = get_first(stack);
		val = first->val;
		if (stack->end != first)
		  stack->end->next = first->next;
		else
		  stack->end = NULL;
		stack_node_destroy(stack->membank, first);
		stack->length -= 1;
	}
	return val;
}

void* stack_exist(const void *const stack_ptr, const void *val,
		  uint32_t (*key)(const void*),
		  int8_t (*compare)(const void*, const void*))
{
	const nb_stack_t *const stack = stack_ptr;
	void *existing_val = NULL;
	if (is_not_empty(stack)) {
		stack_node_t *node = exist_node(stack, val, compare);
		if (NULL != node)
			existing_val = node->val;
	}
	return existing_val;
}

static stack_node_t* exist_node(const nb_stack_t *const stack, const void *val,
			  int8_t (*compare)(const void*, const void*))
{
	stack_node_t *first = get_first(stack);
	stack_node_t *existing_node = NULL;
	if (0 == compare(first->val, val)) {
		existing_node = first;
	} else {
		stack_node_t *iter = first->next;
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
	nb_stack_t *stack = stack_ptr;
	void *deleted_val = NULL;
	if (is_not_empty(stack)) {
		stack_node_t *node = exist_node(stack, val, compare);
		if (NULL != node) {
			unlink_node(stack, node);
			deleted_val = node->val;
			stack_node_destroy(stack->membank, node);
			stack->length -= 1;
		}
	}
	return deleted_val;
}

static void unlink_node(nb_stack_t *stack, const stack_node_t *const node)
{
	stack_node_t *prev = stack_node_get_prev(node);
	if (node == prev /* implies node is equal to stack->end */)
		stack->end = NULL;
	else if (node == stack->end /* implies node is different to prev */)
		stack->end = prev;
	prev->next = node->next;
}

inline uint32_t stack_get_length(const void *const stack_ptr)
{
	const nb_stack_t *const stack = stack_ptr;
	return stack->length;
}

inline bool stack_is_empty(const void *const stack_ptr)
{
	const nb_stack_t *const stack = stack_ptr;
	return (NULL == stack->end);
}

inline bool stack_is_not_empty(const void *const stack_ptr)
{
	return is_not_empty(stack_ptr);
}
