/******************************************************************************
 *   List DST: Circular Linked List                                           *
 *   2011-2015 Victor Eduardo Cardoso Nungaray                                *
 *   Twitter: @victore_cardoso                                                *
 *   email: victorc@cimat.mx                                                  *
 ******************************************************************************/

#include <stdlib.h>
#include <string.h>

#include "queue_node.h"
#include "queue_struct.h"
#include "queue_dst.h"

static bool is_not_empty(const queue_t *const list);
static node_t* get_first(const queue_t *const list);
static void* malloc_queue(void);
static void insert_node_as_starting(queue_t *list, const void *const val);
static void add_node(queue_t *list, node_t *node);
static void add_first_node(queue_t *list, node_t *node);
static void null_destroy(void *val);
static void link_node(queue_t *queue, node_t *node);
static node_t* exist_node(const queue_t *const queue, const void *val,
			  int8_t (*compare)(const void*, const void*));
static void unlink_node(queue_t *queue, const node_t *const node);

inline uint16_t queue_get_memsize(void)
{
	return sizeof(queue_t);
}

inline void queue_init(void *queue_ptr)
{
	queue_t *queue = queue_ptr;
	queue->length = 0;
	queue->end = NULL;
}

void queue_copy(void *queue_ptr, const void *src_queue_ptr,
		void* (*clone)(const void*))
{
	queue_t *queue = queue_ptr;
	const queue_t *src_queue = src_queue_ptr;
	
	if (is_not_empty(src_queue)) {
		queue->length = src_queue->length;
		node_t *end = node_clone(src_queue->end, clone);
		queue->end = end; 

		node_t *iter = get_first(src_queue);
		do {
			node_t *node = node_clone(iter, clone);
			end->next = node;
			end = node;
			iter = iter->next;
		} while (iter != src_queue->end);
		end->next = queue->end;
	}
}


static inline bool is_not_empty(const queue_t *const restrict queue)
{
	return (queue->end != NULL);
}

static inline node_t* get_first(const queue_t *const restrict queue)
{
	return queue->end->next;
}

void queue_clear(void *queue_ptr,
		 void (*destroy)(void*))
{
	queue_t *queue = queue_ptr;
	if (is_not_empty(queue)) {
		node_t* iter = get_first(queue);
		queue->end->next = NULL;
		while (NULL != iter) {
			node_t* to_destroy = iter;
			iter = iter->next;
			node_destroy(to_destroy, destroy);
		}
		queue->length = 0;
		queue->end = NULL;
	}
}

inline void* queue_create(void)
{
	void *queue = malloc_queue();
	queue_init(queue);
	return queue;
}

static inline void* malloc_queue(void)
{
	uint16_t size = queue_get_memsize();
  	return malloc(size);
}

inline void* queue_clone(const void *const queue_ptr,
			 void* (*clone)(const void*))
{
	void *queue = malloc_queue();
	queue_copy(queue, queue_ptr, clone);
	return queue;
}

inline void queue_destroy(void *queue_ptr,
		  void (*destroy)(void*))
{
	queue_clear(queue_ptr, destroy);
	free(queue_ptr);
}

void queue_merge(void *queue1_ptr, void *queue2_ptr,
		 uint32_t (*key)(const void*))
{
	queue_t *queue1 = queue1_ptr;
	queue_t *queue2 = queue2_ptr;
	if (is_not_empty(queue2)) {
		queue1->length += queue2->length;
		if (is_not_empty(queue1)) {
			node_t *first_node = get_first(queue1);
			queue1->end->next = get_first(queue2);
			queue2->end->next = first_node;
		}
		queue1->end = queue2->end;
		queue2->length = 0;
		queue2->end = NULL;
	}
}

inline bool queue_insert_first(void *queue_ptr, const void *val,
			       uint32_t (*key)(const void*))
{
	queue_t *queue = queue_ptr;
	insert_node_as_starting(queue, val);
	return true;
}

static inline void insert_node_as_starting(queue_t *restrict queue,
					   const void *const restrict val)
{
	node_t *const restrict node = node_create();
	node->val = (void*) val;
	add_node(queue, node);
	queue->length += 1;
}

static inline void add_node(queue_t *restrict queue,
			    node_t *restrict node)
{
	if (queue_is_empty(queue))
		add_first_node(queue, node);
	else
		link_node(queue, node);
}

static inline void add_first_node(queue_t *restrict queue,
				  node_t *restrict node)
{
	queue->end = node;
	node->next = node;
}

static inline void link_node(queue_t *restrict queue, 
			     node_t *restrict node)
{
	node->next = get_first(queue);
	queue->end->next = node;
}

inline bool queue_insert(void *queue_ptr, const void *val,
			 uint32_t (*key)(const void*))
{
	queue_t *queue = queue_ptr;
	insert_node_as_starting(queue, val);
	queue->end = queue->end->next;
	return true;
}

inline void* queue_get_first(const void *const queue_ptr)
{
	const queue_t *const queue = queue_ptr;
	void *val = NULL;
	if (is_not_empty(queue)) {
		node_t *first = get_first(queue);
		val = first->val;
	}
	return val;
}

void* queue_delete_first(void *queue_ptr,
			uint32_t (*key)(const void*))
{
	queue_t *queue = queue_ptr;
	void *val = NULL;
	if (is_not_empty(queue)) {
		node_t *first = get_first(queue);
		val = first->val;
		if (queue->end != first)
		  queue->end->next = first->next;
		else
		  queue->end = NULL;
		node_destroy(first, null_destroy);
		queue->length -= 1;
	}
	return val;
}

static inline void null_destroy(void *val)
{
  ; /* Null statement */
}

void* queue_exist(const void *const queue_ptr, const void *val,
		  uint32_t (*key)(const void*),
		  int8_t (*compare)(const void*, const void*))
{
	const queue_t *const queue = queue_ptr;
	void *existing_val = NULL;
	if (is_not_empty(queue)) {
		node_t *node = exist_node(queue, val, compare);
		if (NULL != node)
			existing_val = node->val;
	}
	return existing_val;
}

static node_t* exist_node(const queue_t *const queue, const void *val,
			  int8_t (*compare)(const void*, const void*))
{
	node_t *first = get_first(queue);
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

void* queue_delete(void *queue_ptr, const void *val,
		  uint32_t (*key)(const void*),
		  int8_t (*compare)(const void*, const void*))
{
	queue_t *queue = queue_ptr;
	void *deleted_val = NULL;
	if (is_not_empty(queue)) {
		node_t *node = exist_node(queue, val, compare);
		if (NULL != node) {
			unlink_node(queue, node);
			deleted_val = node->val;
			node_destroy(node, null_destroy);
			queue->length -= 1;
		}
	}
	return deleted_val;
}

static void unlink_node(queue_t *queue, const node_t *const node)
{
	node_t *prev = node_get_prev(node);
	if (node == prev /* implies node is equal to queue->end */)
		queue->end = NULL;
	else if (node == queue->end /* implies node is different to prev */)
		queue->end = prev;
	prev->next = node->next;
}

inline uint32_t queue_get_length(const void *const queue_ptr)
{
	const queue_t *const queue = queue_ptr;
	return queue->length;
}

inline bool queue_is_empty(const void *const queue_ptr)
{
	const queue_t *const queue = queue_ptr;
	return (NULL == queue->end);
}

inline bool queue_is_not_empty(const void *const queue_ptr)
{
	return is_not_empty(queue_ptr);
}
