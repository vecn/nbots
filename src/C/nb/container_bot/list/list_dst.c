/******************************************************************************
 *   List DST: Circular Linked List                                           *
 *   2011-2015 Victor Eduardo Cardoso Nungaray                                *
 *   Twitter: @victore_cardoso                                                *
 *   email: victorc@cimat.mx                                                  *
 ******************************************************************************/

#include <stdlib.h>
#include <string.h>

#include "list_node.h"
#include "list_dst.h"

typedef struct {
	/* Circular linked list */
	uint32_t length;
	node_t *end;
} list_t;

static bool is_not_empty(const list_t *const list);
static node_t* get_first(const list_t *const list);
static void insert_node_as_starting(list_t *list, const void *const val);
static void add_node(list_t *list, node_t *node);
static void add_first_node(list_t *list, node_t *node);
static void null_destroy(void *val);
static void link_node(list_t *list, node_t *node);
static node_t* exist_node(const list_t *const list, const void *const val,
			  bool (*are_equal)(const void *const,
					    const void *const));
static void unlink_node(list_t *list, const node_t *const node);

inline void* list_create(void)
{
  	return calloc(1, sizeof(list_t));
}

void* list_clone(const void *const list_ptr,
		 void* (*clone)(const void *const))
{
	const list_t *const restrict list = list_ptr;
	list_t *const restrict lc = list_create();
	
	if (is_not_empty(list)) {
		lc->length = list->length;
		node_t *end = node_clone(list->end, clone);
		lc->end = end; 

		node_t *iter = get_first(list);
		do {
			node_t *node = node_clone(iter, clone);
			end->next = node;
			end = node;
			iter = iter->next;
		} while (iter != list->end);
		end->next = lc->end;
	}
	return lc;
}

static inline bool is_not_empty(const list_t *const restrict list)
{
	return (list->end != NULL);
}

static inline node_t* get_first(const list_t *const restrict list)
{
	return list->end->next;
}

void list_merge(void *list1_ptr, void *list2_ptr,
		 uint32_t (*key)(const void *const))
{
	list_t *list1 = list1_ptr;
	list_t *list2 = list2_ptr;
	if (is_not_empty(list2)) {
		list1->length += list2->length;
		if (is_not_empty(list1)) {
			node_t *first_node = get_first(list1);
			list1->end->next = get_first(list2);
			list2->end->next = first_node;
		}
		list1->end = list2->end;
		list2->length = 0;
		list2->end = NULL;
	}
}

inline void list_destroy(void *list_ptr,
			 void (*destroy)(void*))
{
	list_clear(list_ptr, destroy);
	free(list_ptr);
}

void list_clear(void *list_ptr,
		 void (*destroy)(void*))
{
	list_t *list = list_ptr;
	if (is_not_empty(list)) {
		node_t* iter = get_first(list);
		list->end->next = NULL;
		while (NULL != iter) {
			node_t* to_destroy = iter;
			iter = iter->next;
			node_destroy(to_destroy, destroy);
		}
		list->length = 0;
		list->end = NULL;
	}
}

inline bool list_insert_first(void *list_ptr, 
			      const void *const restrict val,
			      uint32_t (*key)(const void *const))
{
	list_t *list = list_ptr;
	insert_node_as_starting(list, val);
	return true;
}

static inline void insert_node_as_starting(list_t *restrict list,
					   const void *const restrict val)
{
	node_t *const restrict node = node_create();
	node->val = (void*) val;
	add_node(list, node);
	list->length += 1;
}

static inline void add_node(list_t *restrict list,
			    node_t *restrict node)
{
	if (list_is_empty(list))
		add_first_node(list, node);
	else
		link_node(list, node);
}

static inline void add_first_node(list_t *restrict list,
				  node_t *restrict node)
{
	list->end = node;
	node->next = node;
}

static inline void link_node(list_t *restrict list, 
			     node_t *restrict node)
{
	node->next = get_first(list);
	list->end->next = node;
}

inline bool list_insert_last(void *list_ptr,
			     const void *const restrict val,
			     uint32_t (*key)(const void *const))
{
	list_t *list = list_ptr;
	insert_node_as_starting(list, val);
	list->end = list->end->next;
	return true;
}

inline void* list_get_first(const void *const list_ptr)
{
	const list_t *const list = list_ptr;
	void *val = NULL;
	if (is_not_empty(list)) {
		node_t *first = get_first(list);
		val = first->val;
	}
	return val;
}

void* list_delete_first(void *list_ptr,
			uint32_t (*key)(const void *const))
{
	list_t *list = list_ptr;
	void *val = NULL;
	if (is_not_empty(list)) {
		node_t *first = get_first(list);
		val = first->val;
		if (list->end != first)
		  list->end->next = first->next;
		else
		  list->end = NULL;
		node_destroy(first, null_destroy);
		list->length -= 1;
	}
	return val;
}

static inline void null_destroy(void *val)
{
  ; /* Null statement */
}

void* list_exist(const void *const list_ptr, const void *const val,
		  uint32_t (*key)(const void *const),
		  bool (*are_equal)(const void *const, const void *const))
{
	const list_t *const list = list_ptr;
	void *existing_val = NULL;
	if (is_not_empty(list)) {
		node_t *node = exist_node(list, val, are_equal);
		if (NULL != node)
			existing_val = node->val;
	}
	return existing_val;
}

static node_t* exist_node(const list_t *const list, const void *const val,
			  bool (*are_equal)(const void *const, 
					    const void *const))
{
	node_t *first = get_first(list);
	node_t *existing_node = NULL;
	if (are_equal(first->val, val)) {
		existing_node = first;
	} else {
		node_t *iter = first->next;
		while (iter != first) {
			if (are_equal(iter->val, val)) {
				existing_node = iter;
				break;
			}
			iter = iter->next;
		}
	}
	return existing_node;
}

void* list_delete(void *list_ptr, const void *const val,
		  uint32_t (*key)(const void *const),
		  bool (*are_equal)(const void *const, const void *const))
{
	list_t *list = list_ptr;
	void *deleted_val = NULL;
	if (is_not_empty(list)) {
		node_t *node = exist_node(list, val, are_equal);
		if (NULL != node) {
			unlink_node(list, node);
			deleted_val = node->val;
			node_destroy(node, null_destroy);
			list->length -= 1;
		}
	}
	return deleted_val;
}

static void unlink_node(list_t *list, const node_t *const node)
{
	node_t *prev = node_get_prev(node);
	if (node == prev /* implies node is equal to list->end */)
		list->end = NULL;
	else if (node == list->end /* implies node is different to prev */)
		list->end = prev;
	prev->next = node->next;
}

inline uint32_t list_get_length(const void *const list_ptr)
{
	const list_t *const list = list_ptr;
	return list->length;
}

inline bool list_is_empty(const void *const list_ptr)
{
	const list_t *const list = list_ptr;
	return (NULL == list->end);
}

inline bool list_is_not_empty(const void *const list_ptr)
{
  return is_not_empty(list_ptr);
}

inline const void* list_get_iterator_start(const void *const list_ptr)
{
	const list_t *const list = list_ptr;
	return list->end->next;
}
