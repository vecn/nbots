/******************************************************************************
 *   List DST: Circular Linked List                                           *
 *   2011-2015 Victor Eduardo Cardoso Nungaray                                *
 *   Twitter: @victore_cardoso                                                *
 *   email: victorc@cimat.mx                                                  *
 ******************************************************************************/

#include <stdlib.h>
#include <string.h>

#include "nb/memory_bot.h"

#include "nb_queue_node.h"
#include "nb_queue_struct.h"
#include "nb_queue_dst.h"

static bool is_not_empty(const nb_queue_t *const list);
static nb_queue_node_t* get_first(const nb_queue_t *const list);
static void* calloc_queue(void);
static void destroy_values(nb_queue_t *queue, void (*destroy)(void*));
static void insert_node_as_starting(nb_queue_t *list, const void *const val);
static void add_node(nb_queue_t *list, nb_queue_node_t *node);
static void add_first_node(nb_queue_t *list, nb_queue_node_t *node);
static void link_node(nb_queue_t *queue, nb_queue_node_t *node);
static nb_queue_node_t* exist_node(const nb_queue_t *const queue,
				   const void *val,
				   int8_t (*compare)(const void*, const void*));
static void unlink_node(nb_queue_t *queue, const nb_queue_node_t *const node);

uint32_t nb_queue_get_memsize(void)
{
	return sizeof(nb_queue_t) + nb_membank_get_memsize();
}

void nb_queue_init(void *queue_ptr)
{
	char *memblock = queue_ptr;
	nb_queue_t *queue = (void*) memblock;
	queue->length = 0;
	queue->end = NULL;
	queue->membank = (void*) (memblock + sizeof(nb_queue_t));
	nb_membank_init(queue->membank, nb_queue_node_get_memsize());
}

void nb_queue_copy(void *queue_ptr, const void *src_queue_ptr,
		   void* (*clone)(const void*))
{
	nb_queue_t *queue = queue_ptr;
	const nb_queue_t *src_queue = src_queue_ptr;

	nb_queue_init(queue);

	queue->length = src_queue->length;

	if (is_not_empty(src_queue)) {
		nb_queue_node_t *end = 
			nb_queue_node_clone(queue->membank,
					    src_queue->end, clone);
		queue->end = end; 

		nb_queue_node_t *iter = get_first(src_queue);
		do {
			nb_queue_node_t *node = 
				nb_queue_node_clone(queue->membank,
						    iter, clone);
			end->next = node;
			end = node;
			iter = iter->next;
		} while (iter != src_queue->end);
		end->next = queue->end;
	} else {
		queue->end = NULL;
	}
}


static inline bool is_not_empty(const nb_queue_t *const restrict queue)
{
	return (queue->end != NULL);
}

static inline nb_queue_node_t* get_first(const nb_queue_t *const restrict queue)
{
	return queue->end->next;
}

void nb_queue_finish(void *queue_ptr, void (*destroy)(void*))
{
	nb_queue_t *queue = queue_ptr;
	destroy_values(queue, destroy);
	nb_membank_finish(queue->membank);
}

void* nb_queue_create(void)
{
	void *queue = calloc_queue();
	nb_queue_init(queue);
	return queue;
}

static inline void* calloc_queue(void)
{
	uint32_t size = nb_queue_get_memsize();
  	return nb_calloc(size);
}

void* nb_queue_clone(const void *const queue_ptr,
		     void* (*clone)(const void*))
{
	void *queue = calloc_queue();
	nb_queue_copy(queue, queue_ptr, clone);
	return queue;
}

void nb_queue_destroy(void *queue_ptr,
		      void (*destroy)(void*))
{
	nb_queue_finish(queue_ptr, destroy);
	free(queue_ptr);
}

void nb_queue_clear(void *queue_ptr,
		    void (*destroy)(void*))
{
	nb_queue_t *queue = queue_ptr;
	destroy_values(queue, destroy);
	nb_membank_clear(queue->membank);
	queue->length = 0;
	queue->end = NULL;
}

static void destroy_values(nb_queue_t *queue, void (*destroy)(void*))
{
	if (NULL != destroy) {
		if (is_not_empty(queue)) {
			nb_queue_node_t* iter = get_first(queue);
			queue->end->next = NULL;
			while (NULL != iter) {
				destroy(iter->val);
				iter = iter->next;
			}
		}
	}
}

void nb_queue_merge(void *queue1_ptr, void *queue2_ptr,
		    uint32_t (*key)(const void*),
		    int8_t (*compare)(const void*, const void*))
{
	nb_queue_t *queue1 = queue1_ptr;
	nb_queue_t *queue2 = queue2_ptr;
	if (is_not_empty(queue2)) {
		queue1->length += queue2->length;
		if (is_not_empty(queue1)) {
			nb_queue_node_t *first_node = get_first(queue1);
			queue1->end->next = get_first(queue2);
			queue2->end->next = first_node;
		}
		queue1->end = queue2->end;
		queue2->length = 0;
		queue2->end = NULL;

		nb_membank_merge(queue1->membank, queue2->membank);
	}
}

bool nb_queue_insert_first(void *queue_ptr, const void *val,
			   uint32_t (*key)(const void*),
			   int8_t (*compare)(const void*, const void*))
{
	nb_queue_t *queue = queue_ptr;
	insert_node_as_starting(queue, val);
	return true;
}

static inline void insert_node_as_starting(nb_queue_t *restrict queue,
					   const void *const restrict val)
{
	nb_queue_node_t *const node = nb_queue_node_create(queue->membank);
	node->val = (void*) val;
	add_node(queue, node);
	queue->length += 1;
}

static inline void add_node(nb_queue_t *restrict queue,
			    nb_queue_node_t *restrict node)
{
	if (nb_queue_is_empty(queue))
		add_first_node(queue, node);
	else
		link_node(queue, node);
}

static inline void add_first_node(nb_queue_t *restrict queue,
				  nb_queue_node_t *restrict node)
{
	queue->end = node;
	node->next = node;
}

static inline void link_node(nb_queue_t *restrict queue, 
			     nb_queue_node_t *restrict node)
{
	node->next = get_first(queue);
	queue->end->next = node;
}

bool nb_queue_insert(void *queue_ptr, const void *val,
		     uint32_t (*key)(const void*),
		     int8_t (*compare)(const void*, const void*))
{
	nb_queue_t *queue = queue_ptr;
	insert_node_as_starting(queue, val);
	queue->end = queue->end->next;
	return true;
}

void* nb_queue_get_first(const void *const queue_ptr)
{
	const nb_queue_t *const queue = queue_ptr;
	void *val = NULL;
	if (is_not_empty(queue)) {
		nb_queue_node_t *first = get_first(queue);
		val = first->val;
	}
	return val;
}

void* nb_queue_delete_first(void *queue_ptr,
			    uint32_t (*key)(const void*))
{
	nb_queue_t *queue = queue_ptr;
	void *val = NULL;
	if (is_not_empty(queue)) {
		nb_queue_node_t *first = get_first(queue);
		val = first->val;
		if (queue->end != first)
			queue->end->next = first->next;
		else
			queue->end = NULL;
		nb_queue_node_destroy(queue->membank, first);
		queue->length -= 1;
	}
	return val;
}

void* nb_queue_exist(const void *const queue_ptr, const void *val,
		     uint32_t (*key)(const void*),
		     int8_t (*compare)(const void*, const void*))
{
	const nb_queue_t *const queue = queue_ptr;
	void *existing_val = NULL;
	if (is_not_empty(queue)) {
		nb_queue_node_t *node = exist_node(queue, val, compare);
		if (NULL != node)
			existing_val = node->val;
	}
	return existing_val;
}

static nb_queue_node_t* exist_node(const nb_queue_t *const queue,
				   const void *val,
				   int8_t (*compare)(const void*, const void*))
{
	nb_queue_node_t *first = get_first(queue);
	nb_queue_node_t *existing_node = NULL;
	if (0 == compare(first->val, val)) {
		existing_node = first;
	} else {
		nb_queue_node_t *iter = first->next;
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

void* nb_queue_delete(void *queue_ptr, const void *val,
		      uint32_t (*key)(const void*),
		      int8_t (*compare)(const void*, const void*))
{
	nb_queue_t *queue = queue_ptr;
	void *deleted_val = NULL;
	if (is_not_empty(queue)) {
		nb_queue_node_t *node = exist_node(queue, val, compare);
		if (NULL != node) {
			unlink_node(queue, node);
			deleted_val = node->val;
			nb_queue_node_destroy(queue->membank, node);
			queue->length -= 1;
		}
	}
	return deleted_val;
}

static void unlink_node(nb_queue_t *queue, const nb_queue_node_t *const node)
{
	nb_queue_node_t *prev = nb_queue_node_get_prev(node);
	if (node == prev /* implies node is equal to queue->end */)
		queue->end = NULL;
	else if (node == queue->end /* implies node is different to prev */)
		queue->end = prev;
	prev->next = node->next;
}

uint32_t nb_queue_get_length(const void *const queue_ptr)
{
	const nb_queue_t *const queue = queue_ptr;
	return queue->length;
}

bool nb_queue_is_empty(const void *const queue_ptr)
{
	const nb_queue_t *const queue = queue_ptr;
	return (NULL == queue->end);
}

bool nb_queue_is_not_empty(const void *const queue_ptr)
{
	return is_not_empty(queue_ptr);
}
