/******************************************************************************
 *   Generic Container                                                        *
 *   2011-2015 Victor Eduardo Cardoso Nungaray                                *
 *   Twitter: @victore_cardoso                                                *
 *   email: victorc@cimat.mx                                                  *
 ******************************************************************************/

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>

#include "list/list_dst.h"
#include "avl/avl_dst.h"
#include "htable/htable_dst.h"
#include "heap/heap_dst.h"

#include "vcn/container_bot/container.h"
#include "vcn/container_bot/iterator.h"

typedef struct {
	uint32_t (*key)(const void *const); 
	void (*destroy)(void*);
	bool (*are_equal)(const void *const, const void *const);
	void* (*clone)(const void *const);
} dst_functions;

struct vcn_container_s {
	int id;
	void *dst;
	dst_functions fdst;
	void* (*create)(void);
	void* (*clone)(const void *const,
		       void* (*clone)(const void *const));
	void (*merge)(void*, void*,
		      uint32_t (*key)(const void *const));
	void (*destroy)(void*,
			void (*destroy)(void*));
	void (*clear)(void*,
		      void (*destroy)(void*));
	bool (*insert)(void*, const void *const,
		       uint32_t (*key)(const void *const));
	void* (*get_first)(const void* const);
	void* (*delete_first)(void*,
			      uint32_t (*key)(const void *const));
	void* (*exist)(const void *const, const void *const,
		       uint32_t (*key)(const void *const),
		       bool (*are_equal)(const void *const, const void *const));
	void* (*delete)(void*, const void *const,
			uint32_t (*key)(const void *const),
			bool (*are_equal)(const void *const, const void *const));
	uint32_t (*get_length)(const void *const);
	bool (*is_empty)(const void *const);
	bool (*is_not_empty)(const void *const);
};

static vcn_container_t *create_without_dst(int8_t id);
static void init_dst_functions(vcn_container_t *container);
static void set_functions(vcn_container_t *container, int8_t id);
static uint32_t key_ptr(const void *const restrict ptr);
static void destroy_null(void* ptr);
static bool are_equal_ptr(const void *const p1, const void *const p2);
static void* clone_same_ptr(const void *const ptr);
static void container_set_queue(vcn_container_t *container);
static void container_set_stack(vcn_container_t *container);
static void container_set_avl(vcn_container_t *container);
static void container_set_heap(vcn_container_t *container);
static void container_set_htable(vcn_container_t *container);
static void container_set_null(vcn_container_t *container);

static void* null_create(void);
static void* null_clone(const void *const obj,
			void* (*clone)(const void *const));
static void null_merge(void *obj1, void *obj2,
		       uint32_t (*key)(const void *const));
static void null_clear_destroy(void *obj,
			       void (*destroy)(void*));
static bool null_insert(void *obj1, const void *const obj2,
			uint32_t (*key)(const void *const));
static void* null_get_first(const void *const obj);
static void* null_delete_first(void *obj,
			       uint32_t (*key)(const void *const));
static void* null_exist(const void *const obj1, const void *const obj2,
			uint32_t (*key)(const void *const),
			bool (*are_equal)(const void *const,
						 const void *const));
static void* null_delete(void *obj1, const void *const obj2,
			 uint32_t (*key)(const void *const),
			 bool (*are_equal)(const void *const,
					   const void *const));
static uint32_t null_get_length(const void *const obj);
static bool null_is_empty(const void *const obj);
static bool null_is_not_empty(const void *const obj);

static void copy_dst_functions(vcn_container_t *dest, 
			       const vcn_container_t *const src);
static bool is_the_same_dst(const vcn_container_t *const c1,
			    const vcn_container_t *const c2);
static void insert_2clear_into_main(vcn_container_t *main, 
				    vcn_container_t *to_clear);
static bool casting_is_valid(int8_t id1, int8_t id2);
static bool is_cast_between_QUEUE_and_STACK(int8_t id1, int8_t id2);
static void cast_container(vcn_container_t* container, int8_t new_id);
static void* queue_do(vcn_container_t *container, const char* func,
		      void *data, int8_t *status);
static void* sorted_do(vcn_container_t *container, const char* func,
		      void *data, int8_t *status);
static void* hash_do(vcn_container_t *container, const char* func,
		      void *data, int8_t *status);

vcn_container_t *vcn_container_create(int8_t id)
{
	vcn_container_t *container = create_without_dst(id);
	container->dst = container->create();
	return container;
}

static vcn_container_t *create_without_dst(int8_t id)
{
	vcn_container_t *container = calloc(1, sizeof(*container));
	container->id = id;
	if (id >= VCN_CONTAINER_NULL)
		container->id = VCN_CONTAINER_NULL;
	init_dst_functions(container);
	set_functions(container, id);
	return container;
}

static void init_dst_functions(vcn_container_t *container)
{
	container->fdst.key = key_ptr;
	container->fdst.destroy = destroy_null;
	container->fdst.are_equal = are_equal_ptr;
	container->fdst.clone = clone_same_ptr;
}

static void set_functions(vcn_container_t *container, int8_t id)
{
	switch (id) {
	case VCN_CONTAINER_QUEUE:
		container_set_queue(container);
		break;
	case VCN_CONTAINER_STACK:
		container_set_stack(container);
		break;
	case VCN_CONTAINER_SORTED:
		container_set_avl(container);
		break;
	case VCN_CONTAINER_HEAP:
		container_set_heap(container);
		break;
	case VCN_CONTAINER_HASH:
		container_set_htable(container);
		break;
	default:
		container_set_null(container);
	}
}

static inline uint32_t key_ptr(const void *const restrict ptr)
{
	return (uint32_t)((uintptr_t)ptr);
}

static inline void destroy_null(void* ptr)
{
	; /* Null statement */
}

static inline bool are_equal_ptr(const void *const p1, const void *const p2)
{
	return (p1 == p2);
}

static void* clone_same_ptr(const void *const ptr)
{
	return (void*) ptr;
}

static void container_set_queue(vcn_container_t *container)
{
	container->create = list_create;
	container->clone = list_clone;
	container->merge = list_merge;
	container->destroy = list_destroy;
	container->clear = list_clear;
	container->insert = list_insert_last;
	container->get_first = list_get_first;
	container->delete_first = list_delete_first;
	container->exist = list_exist;
	container->delete = list_delete;
	container->get_length = list_get_length;
	container->is_empty = list_is_empty;
	container->is_not_empty = list_is_not_empty;
}

static void container_set_stack(vcn_container_t *container)
{
	container->create = list_create;
	container->clone = list_clone;
	container->merge = list_merge;
	container->destroy = list_destroy;
	container->clear = list_clear;
	container->insert = list_insert_first;
	container->get_first = list_get_first;
	container->delete_first = list_delete_first;
	container->exist = list_exist;
	container->delete = list_delete;
	container->get_length = list_get_length;
	container->is_empty = list_is_empty;
	container->is_not_empty = list_is_not_empty;
}

static void container_set_avl(vcn_container_t *container)
{
	container->create = avl_create;
	container->clone = avl_clone;
	container->merge = avl_merge;
	container->destroy = avl_destroy;
	container->clear = avl_clear;
	container->insert = avl_insert;
	container->get_first = avl_get_first;
	container->delete_first = avl_delete_first;
	container->exist = avl_exist;
	container->delete = avl_delete;
	container->get_length = avl_get_length;
	container->is_empty = avl_is_empty;
	container->is_not_empty = avl_is_not_empty;
}

static void container_set_heap(vcn_container_t *container)
{
	container->create = heap_create;
 	container->clone = heap_clone;
	container->merge = heap_merge;
	container->destroy = heap_destroy;
	container->clear = heap_clear;
	container->insert = heap_insert;
	container->get_first = heap_get_first;
	container->delete_first = heap_delete_first;
	container->exist = heap_exist;
	container->delete = heap_delete;
	container->get_length = heap_get_length;
	container->is_empty = heap_is_empty;
	container->is_not_empty = heap_is_not_empty;
}

static void container_set_htable(vcn_container_t *container)
{
	container->create = htable_create;
	container->clone = htable_clone;
	container->merge = htable_merge;
	container->destroy = htable_destroy;
	container->clear = htable_clear;
	container->insert = htable_insert;
	container->get_first = htable_get_first;
	container->delete_first = htable_delete_first;
	container->exist = htable_exist;
	container->delete = htable_delete;
	container->get_length = htable_get_length;
	container->is_empty = htable_is_empty;
	container->is_not_empty = htable_is_not_empty;
}

static void container_set_null(vcn_container_t *container)
{
	container->create = null_create;
	container->clone = null_clone;
	container->merge = null_merge;
	container->destroy = null_clear_destroy;
	container->clear = null_clear_destroy;
	container->insert = null_insert;
	container->get_first = null_get_first;
	container->delete_first = null_delete_first;
	container->exist = null_exist;
	container->delete = null_delete;
	container->get_length = null_get_length;
	container->is_empty = null_is_empty;
	container->is_not_empty = null_is_not_empty;
}

static inline void* null_create(void)
{
	return NULL;
}

static inline void* null_clone(const void *const obj,
			       void* (*clone)(const void *const))
{
  return (void*) obj;
}

static inline void null_merge(void *obj1, void *obj2,
			       uint32_t (*key)(const void *const))
{
  	; /* Null statement */
}

static inline void null_clear_destroy(void *obj,
				      void (*destroy)(void*))
{
  	; /* Null statement */
}

static inline bool null_insert(void* obj1, const void *const obj2,
			       uint32_t (*key)(const void *const))
{
  	return false;
}

static inline void* null_get_first(const void *const obj)
{
  	return NULL;
}

static inline void* null_delete_first(void* obj,
				      uint32_t (*key)(const void *const))
{
  	return NULL;
}

static inline void* null_exist(const void *const obj1, const void *const obj2, 
			       uint32_t (*key)(const void *const),
			       bool (*are_equal)(const void *const,
						 const void *const))
{
  	return NULL;
}

static inline void* null_delete(void *obj1, const void *const obj2, 
				uint32_t (*key)(const void *const),
				bool (*are_equal)(const void *const,
						  const void *const))
{
  	return NULL;
}

static uint32_t null_get_length(const void *const obj)
{
  	return 0;
}

static inline bool null_is_empty(const void *const obj)
{
  	return true;
}

static inline bool null_is_not_empty(const void *const obj)
{
  	return false;
}

vcn_container_t *vcn_container_clone(const vcn_container_t *const container)
{
  	vcn_container_t *cnt_clone = vcn_container_create(container->id);
	cnt_clone->dst = container->clone(container->dst, container->fdst.clone);
	copy_dst_functions(cnt_clone, container);
	return cnt_clone;  
}

static void copy_dst_functions(vcn_container_t *dest,
			       const vcn_container_t *const src)
{
  	dest->fdst.key = src->fdst.key;
	dest->fdst.destroy = src->fdst.destroy;
	dest->fdst.are_equal = src->fdst.are_equal;
	dest->fdst.clone = src->fdst.clone;
}

inline void vcn_container_merge(vcn_container_t *main, 
				vcn_container_t *to_clear)
{	
	if (is_the_same_dst(main, to_clear))
		main->merge(main->dst, to_clear->dst, main->fdst.key);
	else
		insert_2clear_into_main(main, to_clear);
}

inline static bool is_the_same_dst(const vcn_container_t *const restrict c1,
				   const vcn_container_t *const restrict c2)
{
  	return (c1->id == c2->id);
}

static void insert_2clear_into_main(vcn_container_t *main, 
				    vcn_container_t *to_clear)
{
  	while (to_clear->is_not_empty(to_clear->dst)) {
    		void *val = to_clear->delete_first(to_clear->dst, to_clear->fdst.key);
		main->insert(main->dst, val, main->fdst.key);
	}
}

void vcn_container_cast(vcn_container_t* container, int8_t new_id)
{
	if (casting_is_valid(container->id, new_id)) {
		if (is_cast_between_QUEUE_and_STACK(container->id, new_id))
			set_functions(container, new_id);
		else
			cast_container(container, new_id);
		container->id = new_id;
	}
}

static inline bool casting_is_valid(int8_t id1, int8_t id2)
{
	return id1 != id2 && 
		id1 < VCN_CONTAINER_NULL &&
		id2 < VCN_CONTAINER_NULL;
}

static inline bool is_cast_between_QUEUE_and_STACK(int8_t id1, int8_t id2)
{
	return (id1 == VCN_CONTAINER_QUEUE && id2 == VCN_CONTAINER_STACK) ||
		(id1 == VCN_CONTAINER_STACK && id2 == VCN_CONTAINER_QUEUE);
}

static void cast_container(vcn_container_t* container, int8_t new_id)
{
	vcn_container_t container_old;
	set_functions(&container_old, container->id);
	container_old.dst = container->dst;
	set_functions(container, new_id);
	container->dst = container->create();
	while (container_old.is_not_empty(container_old.dst)) {
		void *val = container_old.delete_first(container_old.dst,
						       container->fdst.key);
		container->insert(container->dst, val, container->fdst.key);
	}
	container_old.destroy(container_old.dst, destroy_null);
}

void** vcn_container_cast_to_array(vcn_container_t *container)
{
  	uint32_t N = container->get_length(container->dst);
	void **array = malloc(N * sizeof(*array));

	N = 0;
	while (container->is_not_empty(container->dst)) {
	  	void *val = container->delete_first(container->dst,
						    container->fdst.key);
		array[N++] = val;
	}
	container->destroy(container->dst, destroy_null);
	free(container);

	return array;
}

inline void vcn_container_destroy(vcn_container_t *container)
{
  	container->destroy(container->dst, container->fdst.destroy);
	free(container);
}

inline void vcn_container_clear(vcn_container_t *container)
{
  	container->clear(container->dst, container->fdst.destroy);
}

void vcn_container_copy_to_array(const vcn_container_t *const cont_src,
				 void **array_dest)
{
	vcn_iterator_t *iter = vcn_iterator_create();
	vcn_iterator_set_container(iter, cont_src);
	uint32_t i = 0;
	while (vcn_iterator_has_more(iter)) {
	  void *val = (void*) vcn_iterator_get_next(iter);
	  array_dest[i] = val;
	  i += 1;
	}
	vcn_iterator_destroy(iter);
}

inline void vcn_container_set_key_generator(vcn_container_t *container,
					    uint32_t (*key)(const void *const))
{
  	container->fdst.key = key;
}

inline void vcn_container_set_destroyer(vcn_container_t *container,
					void (*destroy)(void*))
{
  	container->fdst.destroy = destroy;
}

inline void vcn_container_set_comparer(vcn_container_t *container,
				       bool (*are_equal)(const void *const, 
							 const void *const))
{
  	container->fdst.are_equal = are_equal;
}

inline void vcn_container_set_cloner(vcn_container_t *container,
				     void* (*clone)(const void *const))
{
  	container->fdst.clone = clone;
}

inline bool vcn_container_insert(vcn_container_t *container, 
				 const void *const val)
{
	return container->insert(container->dst, val, container->fdst.key);
}

void vcn_container_insert_array(vcn_container_t *container,
				uint32_t N, void **array)
{
  	for (uint32_t i = 0; i < N; i++)
	 	container->insert(container->dst, array[i],
				  container->fdst.key);
}

inline void* vcn_container_get_first(const vcn_container_t *const container)
{
  	return container->get_first(container->dst);
}

inline void* vcn_container_delete_first(vcn_container_t *container)
{
  	return container->delete_first(container->dst, container->fdst.key);
}

inline void* vcn_container_exist(const vcn_container_t *const container,
				 const void *const val)
{
  	return container->exist(container->dst, val, container->fdst.key,
				container->fdst.are_equal);
}

inline void* vcn_container_delete(vcn_container_t *container,
				  const void *const val)
{
  	return container->delete(container->dst, val,
				 container->fdst.key,
				 container->fdst.are_equal);
}

inline uint32_t vcn_container_get_length(const vcn_container_t *const container)
{
  	return container->get_length(container->dst);
}

inline bool vcn_container_is_empty(const vcn_container_t *const container)
{
  	return container->is_empty(container->dst);
}

inline bool vcn_container_is_not_empty(const vcn_container_t *const container)
{
  	return container->is_not_empty(container->dst);
}
  
inline int8_t vcn_container_get_id(const vcn_container_t *const container)
{
  	return container->id;
}

void* vcn_container_do(vcn_container_t *container, const char* func,
		       void *data, int8_t *status)
{
	void *out = NULL;
	switch (container->id) {
	case VCN_CONTAINER_QUEUE:
		out = queue_do(container, func, data, status);
		break;
	case VCN_CONTAINER_SORTED:
		out = sorted_do(container, func, data, status);
		break;
	case VCN_CONTAINER_HASH:
		out = hash_do(container, func, data, status);
		break;
	default:
		*status = 1;
	}
	return out;
}

static void* queue_do(vcn_container_t *container, const char* func,
		      void *data, int8_t *status)
{
	*status = 0;
	void *out = NULL;
	if (0 == strcmp("insert_first", func)) {
		list_insert_first(container->dst, data, container->fdst.key);
	} else {
		*status = 2;
	}
	return out;
}

static void* sorted_do(vcn_container_t *container, const char* func,
		      void *data, int8_t *status)
{
	*status = 0;
	void *out = NULL;
	if (0 == strcmp("delete_last", func)) {
		out = avl_delete_last(container->dst, container->fdst.key);
	} else {
		*status = 2;
	}
	return out;
}

static void* hash_do(vcn_container_t *container, const char* func,
		     void *data, int8_t *status)
{
	*status = 0;
	void *out = NULL;
	if (0 == strcmp("get_size", func)) {
		uint32_t *size = malloc(sizeof(*size));
		*size = htable_get_size(container->dst);
		out = size;
	} else if (0 == strcmp("get_N_collisions", func)) {
		uint32_t *N = malloc(sizeof(*N));
		*N = htable_get_N_collisions(container->dst);
		out = N;
	} else if (0 == strcmp("get_collisions", func)) {
		out = htable_get_collisions(container->dst);
	} else {
		*status = 2;
	}
	return out;
}

inline void* vcn_container_get_dst(const vcn_container_t *const container)
{
  	return container->dst;
}
