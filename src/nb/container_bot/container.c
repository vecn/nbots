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

#include "queue/nb_queue_binding.h"
#include "stack/binding.h"
#include "sorted/binding.h"
#include "heap/binding.h"
#include "hash/binding.h"

#include "nb/container_bot/container.h"
#include "nb/container_bot/iterator.h"

#include "container_struct.h"

static uint16_t dst_get_memsize(nb_container_type type);
static void init_val_operators(nb_container_t *container);
static uint32_t key_ptr(const void *ptr);
static int8_t compare_ptr(const void *p1, const void *p2);
static void* clone_same_ptr(const void *ptr);
static void set_functions(nb_container_t *container,
			  nb_container_type type);
static void copy_val_operators(nb_container_t *dest,
			       const nb_container_t *const src);
static void *malloc_container(nb_container_type type);

static bool is_the_same_dst(const nb_container_t *const c1,
			    const nb_container_t *const c2);
static void insert_2clear_into_main(nb_container_t *main, 
				    nb_container_t *to_clear);

uint16_t nb_container_get_memsize(nb_container_type type)
{
	uint16_t dst_size = dst_get_memsize(type);
	return dst_size + sizeof(nb_container_t);
}

static uint16_t dst_get_memsize(nb_container_type type)
{
	uint16_t dst_size;
	switch(type) {
	case NB_QUEUE:
		dst_size = nb_queue_get_memsize();
		break;
	case NB_STACK:
		dst_size = stack_get_memsize();
		break;
	case NB_SORTED:
		dst_size = avl_get_memsize();
		break;
	case NB_HEAP:
		dst_size = heap_get_memsize();
		break;
	case NB_HASH:
		dst_size = hash_get_memsize();
		break;
	default:
		dst_size = nb_queue_get_memsize();	
	}
	return dst_size;
}

void nb_container_init(void *container_ptr, nb_container_type type)
{

	nb_container_t *container = container_ptr;
	memset(container, 0, nb_container_get_memsize(type));

	container->type = type;
	if (type >= NB_NULL) {
		container->type = NB_QUEUE;
	}
	uint16_t size = sizeof(nb_container_t);
	char *memblock = container_ptr;
	container->cnt = memblock + size;

	init_val_operators(container);
	set_functions(container, type);
	container->b.init(container->cnt);
}


static void init_val_operators(nb_container_t *container)
{
	container->op.key = key_ptr;
	container->op.compare = compare_ptr;
	container->op.clone = clone_same_ptr;
}


static inline uint32_t key_ptr(const void *ptr)
{
	return (uint32_t)((uintptr_t)ptr);
}

static int8_t compare_ptr(const void *p1, const void *p2)
{
	int8_t out;
	if (p1 < p2)
		out = -1;
	else if (p1 > p2)
		out = 1;
	else
		out = 0;
	return out;
}

static void* clone_same_ptr(const void *ptr)
{
	return (void*) ptr;
}

static void set_functions(nb_container_t *container,
			  nb_container_type type)
{
	switch (type) {
	case NB_QUEUE:
		nb_queue_set_handlers(container);
		break;
	case NB_STACK:
		stack_set_handlers(container);
		break;
	case NB_SORTED:
		sorted_set_handlers(container);
		break;
	case NB_HEAP:
		heap_set_handlers(container);
		break;
	case NB_HASH:
		hash_set_handlers(container);
		break;
	default:
		nb_queue_set_handlers(container);
	}
}

void nb_container_copy(void *container_ptr, const void *src_container_ptr)
{
	nb_container_t *container = container_ptr;
	const nb_container_t *src_container = src_container_ptr;
	
	container->type = src_container->type;

	uint16_t size = sizeof(nb_container_t);
	char *memblock = container_ptr;
	container->cnt = memblock + size;

	copy_val_operators(container, src_container);
	set_functions(container, src_container->type);
	container->b.copy(container->cnt, src_container->cnt,
			  src_container->op.clone);
}

static void copy_val_operators(nb_container_t *dest,
			       const nb_container_t *const src)
{
  	dest->op.key = src->op.key;
	dest->op.destroy = src->op.destroy;
	dest->op.compare = src->op.compare;
	dest->op.clone = src->op.clone;
}

inline void nb_container_finish(void* container_ptr)
{
	nb_container_t *container = container_ptr;
	container->b.finish(container->cnt, container->op.destroy);
}

inline void* nb_container_create(nb_container_type type)
{
	void *container = malloc_container(type);
	nb_container_init(container, type);
	return container;
}

static inline void *malloc_container(nb_container_type type)
{
	uint16_t size = nb_container_get_memsize(type);
	return malloc(size);
}

inline void* nb_container_clone(const void *container_ptr)
{
	nb_container_type type = nb_container_get_type(container_ptr);
	void *container = malloc_container(type);
	nb_container_copy(container, container_ptr);
	return container;
}

inline void nb_container_destroy(void *container_ptr)
{
	nb_container_finish(container_ptr);
	free(container_ptr);
}

inline void nb_container_clear(void* container_ptr)
{
	nb_container_t *container = container_ptr;
	container->b.clear(container->cnt, container->op.destroy);
}

inline void nb_container_merge(nb_container_t *main, 
			       nb_container_t *to_clear)
{	
	if (is_the_same_dst(main, to_clear))
		main->b.merge(main->cnt, to_clear->cnt,
			      main->op.key, main->op.compare);
	else
		insert_2clear_into_main(main, to_clear);
}

inline static bool is_the_same_dst(const nb_container_t *const restrict c1,
				   const nb_container_t *const restrict c2)
{
  	return (c1->type == c2->type);
}

static void insert_2clear_into_main(nb_container_t *main, 
				    nb_container_t *to_clear)
{
  	while (to_clear->b.is_not_empty(to_clear->cnt)) {
    		void *val = to_clear->b.delete_first(to_clear->cnt, 
						     to_clear->op.key);
		main->b.insert(main->cnt, val, main->op.key,
			       main->op.compare);
	}
}


void** nb_container_cast_to_array(nb_container_t *container)
{
  	uint32_t N = container->b.get_length(container->cnt);
	void **array = malloc(N * sizeof(*array));

	N = 0;
	while (container->b.is_not_empty(container->cnt)) {
	  	void *val = container->b.delete_first(container->cnt,
						      container->op.key);
		array[N++] = val;
	}
	nb_container_destroy(container);
	return array;
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

void nb_container_set_key_generator(nb_container_t *container,
				    uint32_t (*key)(const void*))
{
  	container->op.key = key;
}

void nb_container_set_destroyer(nb_container_t *container,
				void (*destroy)(void*))
{
  	container->op.destroy = destroy;
}

void nb_container_unset_destroyer(nb_container_t *container)
{
	container->op.destroy = NULL;
}

void nb_container_set_comparer(nb_container_t *container,
			       int8_t (*compare)(const void*, 
						 const void*))
{
  	container->op.compare = compare;
}

void nb_container_set_cloner(nb_container_t *container,
			     void* (*clone)(const void*))
{
  	container->op.clone = clone;
}

bool nb_container_insert(nb_container_t *container, const void *val)
{
	return container->b.insert(container->cnt, val, container->op.key,
				   container->op.compare);
}

void nb_container_insert_array(nb_container_t *container,
				uint32_t N, void **array)
{
  	for (uint32_t i = 0; i < N; i++)
	 	container->b.insert(container->cnt, array[i],
				    container->op.key,
				    container->op.compare);
}

inline void* nb_container_get_first(const nb_container_t *const container)
{
  	return container->b.get_first(container->cnt);
}

inline void* nb_container_delete_first(nb_container_t *container)
{
  	return container->b.delete_first(container->cnt, container->op.key);
}

inline void* nb_container_exist(const nb_container_t *const container,
				 const void *const val)
{
  	return container->b.exist(container->cnt, val, container->op.key,
				  container->op.compare);
}

inline void* nb_container_delete(nb_container_t *container,
				  const void *const val)
{
  	return container->b.delete(container->cnt, val,
				   container->op.key,
				   container->op.compare);
}

inline uint32_t nb_container_get_length(const nb_container_t *const container)
{
  	return container->b.get_length(container->cnt);
}

inline bool nb_container_is_empty(const nb_container_t *const container)
{
  	return container->b.is_empty(container->cnt);
}

inline bool nb_container_is_not_empty(const nb_container_t *const container)
{
  	return container->b.is_not_empty(container->cnt);
}
  
inline nb_container_type nb_container_get_type
				(const nb_container_t *const container)
{
  	return container->type;
}

void* nb_container_do(nb_container_t *container, const char* func,
		       void *data, int8_t *status)
{
	void *out = NULL;
	switch (container->type) {
	case NB_QUEUE:
		out = nb_queue_do(container, func, data, status);
		break;
	case NB_STACK:
		out = stack_do(container, func, data, status);
		break;
	case NB_SORTED:
		out = sorted_do(container, func, data, status);
		break;
	case NB_HEAP:
		out = heap_do(container, func, data, status);
		break;
	case NB_HASH:
		out = hash_do(container, func, data, status);
		break;
	default:
		*status = 1;
	}
	return out;
}
