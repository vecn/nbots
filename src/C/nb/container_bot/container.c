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

#include "nb/container_bot/container.h"
#include "nb/container_bot/iterator.h"

typedef struct {
	uint32_t (*key)(const void*); 
	void (*destroy)(void*);
	int8_t (*compare)(const void*, const void*);
	void* (*clone)(const void*);
} dst_functions;

struct nb_container_s {
	int id;
	void *dst;
	dst_functions fdst;
	uint16_t (*get_memsize)(void);
	void (*init)(void*);
	void (*copy)(void*, const void*);
	void (*clear)(void*,
		      void (*destroy)(void*)|);
	void* (*create)(void);
	void* (*clone)(const void *,
		       void* (*clone)(const void*));
	void (*destroy)(void*,
			void (*destroy)(void*));
	void (*merge)(void*, void*,
		      uint32_t (*key)(const void*));
	bool (*insert)(void*, const void *const,
		       uint32_t (*key)(const void*));
	void* (*get_first)(const void* const);
	void* (*delete_first)(void*,
			      uint32_t (*key)(const void*));
	void* (*exist)(const void*, const void*,
		       uint32_t (*key)(const void*),
		       int8_t (*compare)(const void*, const void*));
	void* (*delete)(void*, const void*,
			uint32_t (*key)(const void*),
			int8_t (*compare)(const void*, const void*));
	uint32_t (*get_length)(const void*);
	bool (*is_empty)(const void*);
	bool (*is_not_empty)(const void*);
};

static nb_container_t *create_without_dst(int8_t id);
static void init_dst_functions(nb_container_t *container);
static void set_functions(nb_container_t *container, int8_t id);
static uint32_t key_ptr(const void *const restrict ptr);
static void destroy_null(void* ptr);
static bool are_equal_ptr(const void *const p1, const void *const p2);
static void* clone_same_ptr(const void *const ptr);
static void container_set_queue(nb_container_t *container);
static void container_set_stack(nb_container_t *container);
static void container_set_avl(nb_container_t *container);
static void container_set_heap(nb_container_t *container);
static void container_set_htable(nb_container_t *container);
static void container_set_null(nb_container_t *container);

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

static void copy_dst_functions(nb_container_t *dest, 
			       const nb_container_t *const src);
static bool is_the_same_dst(const nb_container_t *const c1,
			    const nb_container_t *const c2);
static void insert_2clear_into_main(nb_container_t *main, 
				    nb_container_t *to_clear);
static bool casting_is_valid(int8_t id1, int8_t id2);
static bool is_cast_between_QUEUE_and_STACK(int8_t id1, int8_t id2);
static void cast_container(nb_container_t* container, int8_t new_id);
static void* queue_do(nb_container_t *container, const char* func,
		      void *data, int8_t *status);
static void* sorted_do(nb_container_t *container, const char* func,
		      void *data, int8_t *status);
static void* hash_do(nb_container_t *container, const char* func,
		      void *data, int8_t *status);

uint16_t nb_container_get_memsize(nb_container_type type)
{
	uint16_t dst_size;
	switch(type) {
	case NB_QUEUE:
		dst_size = list_get_memsize();
		break;
	case NB_STACK:
		dst_size = list_get_memsize();
		break;
	case NB_SORTED:
		dst_size = avl_get_memsize();
		break;
	case NB_HEAP:
		dst_size = heap_get_memsize();
		break;
	case NB_HASH:
		dst_size = htable_get_memsize();
		break;
	default:
		dst_size = list_get_memsize();	
	}
	return dst_size + sizeof(nb_container_t);
}

void nb_container_init(void *container_ptr, nb_container_type type);
void nb_container_copy(void *container_ptr, const void *src_container_ptr);
void nb_container_clear(void* container_ptr);
void* nb_container_create(nb_container_type type);
void* nb_container_clone(const void *container_ptr);
void nb_container_destroy(void *container_ptr);

nb_container_t *nb_container_create(int8_t id)
{
	nb_container_t *container = create_without_dst(id);
	container->dst = container->create();
	return container;
}

static nb_container_t *create_without_dst(int8_t id)
{
	nb_container_t *container = calloc(1, sizeof(*container));
	container->id = id;
	if (id >= NB_NULL)
		container->id = NB_NULL;
	init_dst_functions(container);
	set_functions(container, id);
	return container;
}

static void init_dst_functions(nb_container_t *container)
{
	container->fdst.key = key_ptr;
	container->fdst.destroy = destroy_null;
	container->fdst.are_equal = are_equal_ptr;
	container->fdst.clone = clone_same_ptr;
}

static void set_functions(nb_container_t *container, int8_t id)
{
	switch (id) {
	case NB_QUEUE:
		container_set_queue(container);
		break;
	case NB_STACK:
		container_set_stack(container);
		break;
	case NB_SORTED:
		container_set_avl(container);
		break;
	case NB_HEAP:
		container_set_heap(container);
		break;
	case NB_HASH:
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

static void container_set_queue(nb_container_t *container)
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

static void container_set_stack(nb_container_t *container)
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

static void container_set_avl(nb_container_t *container)
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

static void container_set_heap(nb_container_t *container)
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

static void container_set_htable(nb_container_t *container)
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

static void container_set_null(nb_container_t *container)
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

nb_container_t *nb_container_clone(const nb_container_t *const container)
{
  	nb_container_t *cnt_clone = nb_container_create(container->id);
	cnt_clone->dst = container->clone(container->dst, container->fdst.clone);
	copy_dst_functions(cnt_clone, container);
	return cnt_clone;  
}

static void copy_dst_functions(nb_container_t *dest,
			       const nb_container_t *const src)
{
  	dest->fdst.key = src->fdst.key;
	dest->fdst.destroy = src->fdst.destroy;
	dest->fdst.are_equal = src->fdst.are_equal;
	dest->fdst.clone = src->fdst.clone;
}

inline void nb_container_merge(nb_container_t *main, 
				nb_container_t *to_clear)
{	
	if (is_the_same_dst(main, to_clear))
		main->merge(main->dst, to_clear->dst, main->fdst.key);
	else
		insert_2clear_into_main(main, to_clear);
}

inline static bool is_the_same_dst(const nb_container_t *const restrict c1,
				   const nb_container_t *const restrict c2)
{
  	return (c1->id == c2->id);
}

static void insert_2clear_into_main(nb_container_t *main, 
				    nb_container_t *to_clear)
{
  	while (to_clear->is_not_empty(to_clear->dst)) {
    		void *val = to_clear->delete_first(to_clear->dst, to_clear->fdst.key);
		main->insert(main->dst, val, main->fdst.key);
	}
}

void nb_container_cast(nb_container_t* container, int8_t new_id)
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
		id1 < NB_NULL &&
		id2 < NB_NULL;
}

static inline bool is_cast_between_QUEUE_and_STACK(int8_t id1, int8_t id2)
{
	return (id1 == NB_QUEUE && id2 == NB_STACK) ||
		(id1 == NB_STACK && id2 == NB_QUEUE);
}

static void cast_container(nb_container_t* container, int8_t new_id)
{
	nb_container_t container_old;
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

void** nb_container_cast_to_array(nb_container_t *container)
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

inline void nb_container_destroy(nb_container_t *container)
{
  	container->destroy(container->dst, container->fdst.destroy);
	free(container);
}

inline void nb_container_clear(nb_container_t *container)
{
  	container->clear(container->dst, container->fdst.destroy);
}

void nb_container_copy_to_array(const nb_container_t *const cont_src,
				 void **array_dest)
{
	nb_iterator_t *iter = nb_iterator_create();
	nb_iterator_set_container(iter, cont_src);
	uint32_t i = 0;
	while (nb_iterator_has_more(iter)) {
	  void *val = (void*) nb_iterator_get_next(iter);
	  array_dest[i] = val;
	  i += 1;
	}
	nb_iterator_destroy(iter);
}

inline void nb_container_set_key_generator(nb_container_t *container,
					    uint32_t (*key)(const void *const))
{
  	container->fdst.key = key;
}

inline void nb_container_set_destroyer(nb_container_t *container,
					void (*destroy)(void*))
{
  	container->fdst.destroy = destroy;
}

inline void nb_container_set_comparer(nb_container_t *container,
				       bool (*are_equal)(const void *const, 
							 const void *const))
{
  	container->fdst.are_equal = are_equal;
}

inline void nb_container_set_cloner(nb_container_t *container,
				     void* (*clone)(const void *const))
{
  	container->fdst.clone = clone;
}

inline bool nb_container_insert(nb_container_t *container, 
				 const void *const val)
{
	return container->insert(container->dst, val, container->fdst.key);
}

void nb_container_insert_array(nb_container_t *container,
				uint32_t N, void **array)
{
  	for (uint32_t i = 0; i < N; i++)
	 	container->insert(container->dst, array[i],
				  container->fdst.key);
}

inline void* nb_container_get_first(const nb_container_t *const container)
{
  	return container->get_first(container->dst);
}

inline void* nb_container_delete_first(nb_container_t *container)
{
  	return container->delete_first(container->dst, container->fdst.key);
}

inline void* nb_container_exist(const nb_container_t *const container,
				 const void *const val)
{
  	return container->exist(container->dst, val, container->fdst.key,
				container->fdst.are_equal);
}

inline void* nb_container_delete(nb_container_t *container,
				  const void *const val)
{
  	return container->delete(container->dst, val,
				 container->fdst.key,
				 container->fdst.are_equal);
}

inline uint32_t nb_container_get_length(const nb_container_t *const container)
{
  	return container->get_length(container->dst);
}

inline bool nb_container_is_empty(const nb_container_t *const container)
{
  	return container->is_empty(container->dst);
}

inline bool nb_container_is_not_empty(const nb_container_t *const container)
{
  	return container->is_not_empty(container->dst);
}
  
inline int8_t nb_container_get_id(const nb_container_t *const container)
{
  	return container->id;
}

void* nb_container_do(nb_container_t *container, const char* func,
		       void *data, int8_t *status)
{
	void *out = NULL;
	switch (container->id) {
	case NB_QUEUE:
		out = queue_do(container, func, data, status);
		break;
	case NB_SORTED:
		out = sorted_do(container, func, data, status);
		break;
	case NB_HASH:
		out = hash_do(container, func, data, status);
		break;
	default:
		*status = 1;
	}
	return out;
}

static void* queue_do(nb_container_t *container, const char* func,
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

static void* sorted_do(nb_container_t *container, const char* func,
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

static void* hash_do(nb_container_t *container, const char* func,
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

inline void* nb_container_get_dst(const nb_container_t *const container)
{
  	return container->dst;
}
