/******************************************************************************
 *   Iterator                                                                 *
 *   2011-2015 Victor Eduardo Cardoso Nungaray                                *
 *   Twitter: @victore_cardoso                                                *
 *   email: victorc@cimat.mx                                                  *
 ******************************************************************************/

#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>

#include "queue/binding.h"
#include "stack/binding.h"
#include "sorted/binding.h"
#include "heap/binding.h"
#include "hash/binding.h"

#include "nb/container_bot/container.h"
#include "nb/container_bot/iterator.h"

#include "container_struct.h"
#include "iterator_struct.h"

#define MAX(a,b) (((a)>(b))?(a):(b))

static uint16_t dst_iter_get_max_memsize(void);
static void null_init(void *iter);
static void null_copy(void *iter, const void *src_iter);
static void null_clear(void *iter);
static void* null_create(void);
static void* null_clone(const void *iter);
static void null_destroy(void *iter);
static void null_restart(void *iter);
static const void* null_get_next(void *iter);
static bool null_has_more(const void *const iter);
static void null_set_dst(void *iter, const void *dst);
static void copy_handlers(nb_iterator_t *iter, const nb_iterator_t *src_iter);
static void* malloc_iterator(void);

uint16_t nb_iterator_get_memsize(void)
{
	uint16_t size = dst_iter_get_max_memsize();
	return size + sizeof(nb_iterator_t);
}

static uint16_t dst_iter_get_max_memsize(void)
{
	uint16_t size = MAX(queue_iter_get_memsize(),
			    stack_iter_get_memsize());
	size = MAX(size, avl_iter_get_memsize());
	size = MAX(size, heap_iter_get_memsize());
	size = MAX(size, hash_iter_get_memsize());
	return size;
}

void nb_iterator_init(void *iter_ptr)
{
	nb_iterator_t *iter = iter_ptr;
	iter->b.init = null_init;
	iter->b.copy = null_copy;
	iter->b.clear = null_clear;
	iter->b.create = null_create;
	iter->b.clone = null_clone;
	iter->b.destroy = null_destroy;
	iter->b.restart = null_restart;
	iter->b.get_next = null_get_next;
	iter->b.has_more = null_has_more;
	iter->b.set_dst = null_set_dst;

	char *memblock = iter_ptr;
	iter->internal = memblock + sizeof(nb_iterator_t);
}

static void null_init(void *iter)
{
	; /* Null statement */
}

static void null_copy(void *iter, const void *src_iter)
{
	; /* Null statement */
}

static void null_clear(void *iter)
{
	; /* Null statement */
}

static void* null_create(void)
{
	return NULL;
}

static inline void* null_clone(const void *const iter)
{
	return NULL;
}

static inline void null_destroy(void *iter)
{
	; /* Null statement */
}

static inline void null_restart(void *iter)
{
	; /*  Null statement */
}

static inline const void* null_get_next(void *iter)
{
	return NULL;
}

static inline bool null_has_more(const void *const iter)
{
	return false;
}

static void null_set_dst(void *iter, const void *dst)
{
	; /* Null statement */
}

void nb_iterator_copy(void *iter_ptr, const void *src_iter_ptr)
{
	nb_iterator_t *iter = iter_ptr;
	const nb_iterator_t *src_iter = src_iter_ptr;
	copy_handlers(iter, src_iter);
	iter->b.copy(iter->internal, src_iter->internal);
}

static void copy_handlers(nb_iterator_t *iter, const nb_iterator_t *src_iter)
{
	iter->b.init = src_iter->b.init;
	iter->b.copy = src_iter->b.copy;
	iter->b.clear = src_iter->b.clear;
	iter->b.create = src_iter->b.create;
	iter->b.clone = src_iter->b.clone;
	iter->b.destroy = src_iter->b.destroy;
	iter->b.restart = src_iter->b.restart;
	iter->b.get_next = src_iter->b.get_next;
	iter->b.has_more = src_iter->b.has_more;
}

inline void nb_iterator_clear(void* iter_ptr)
{
	nb_iterator_t *iter = iter_ptr;
	if (NULL != iter->internal)
		iter->b.clear(iter->internal);
}

inline void* nb_iterator_create(void)
{
	void *iterator = malloc_iterator();
	nb_iterator_init(iterator);
	return iterator;
}

static inline void* malloc_iterator(void)
{
	uint16_t size = nb_iterator_get_memsize();
	return malloc(size);
}

inline void* nb_iterator_clone(const void *iter_ptr)
{
	void *iterator = malloc_iterator();
	nb_iterator_copy(iterator, iter_ptr);
	return iterator;
}

inline void nb_iterator_destroy(void *iter_ptr)
{
	nb_iterator_clear(iter_ptr);
	free(iter_ptr);
}



void nb_iterator_set_container(nb_iterator_t *iter,
				const nb_container_t *const container)
{
	switch (container->type) {
	case NB_QUEUE:
		queue_iterator_set_handlers(iter);
		break;
	case NB_STACK:
		stack_iterator_set_handlers(iter);
		break;
	case NB_SORTED:
		sorted_iterator_set_handlers(iter);
		break;
	case NB_HEAP:
		heap_iterator_set_handlers(iter);
		break;
	case NB_HASH:
		hash_iterator_set_handlers(iter);
		break;
	case NB_NULL:
		; /* Null statement */
	}
	iter->b.init(iter->internal);
	iter->b.set_dst(iter->internal, container->cnt);
}

inline void nb_iterator_restart(nb_iterator_t *iter)
{
	iter->b.restart(iter->internal);
}

inline const void* nb_iterator_get_next(nb_iterator_t *iter)
{
	return iter->b.get_next(iter->internal);
}

inline bool nb_iterator_has_more(const nb_iterator_t *const iter)
{
	return iter->b.has_more(iter->internal);
}
