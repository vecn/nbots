/******************************************************************************
 *   Iterator                                                                 *
 *   2011-2015 Victor Eduardo Cardoso Nungaray                                *
 *   Twitter: @victore_cardoso                                                *
 *   email: victorc@cimat.mx                                                  *
 ******************************************************************************/

#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>

#include "list/list_iterator.h"
#include "avl/avl_iterator.h"
#include "htable/htable_iterator.h"
#include "heap/heap_iterator.h"

#include "nb/container_bot/container.h"
#include "nb/container_bot/iterator.h"

#include "container_struct.h"

#define MAX(a,b) (((a)>(b))?(a):(b))

struct nb_iterator_s {
	void *dst_iter;
	uint16_t (*get_memsize)(void);
	void (*init)(void*);
	void (*copy)(void*, const void*);
	void (*clear)(void*);
	void* (*create)(void);
	void* (*clone)(const void *);
	void (*destroy)(void*);
	void (*restart)(void*);
	const void* (*get_next)(void*);
	bool (*has_more)(const void *);
	void (*set_dst)(void*, const void*)
};

static uint16_t dst_iter_get_max_memsize(void);
static uint16_t null_get_memsize(void);
static void null_init(void *iter);
static void null_copy(void *iter, const void *src_iter);
static void null_clear(void *iter);
static void null_create(void);
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

static uint16_t dst_iter_get_memsize(void)
{
	uint16_t size = MAX(queue_iter_get_memsize(),
			    stack_iter_get_memsize());
	size = MAX(size, sorted_iter_get_memsize());
	size = MAX(size, heap_iter_get_memsize());
	size = MAX(size, hash_iter_get_memsize());
	return size;
}

void nb_iterator_init(void *iter_ptr)
{
	nb_iterator_t *iter = iter_ptr;
	iter->get_memsize = null_get_memsize;
	iter->init = null_init;
	iter->copy = null_copy;
	iter->clear = null_clear;
	iter->create = null_create;
	iter->clone = null_clone;
	iter->destroy = null_destroy;
	iter->restart = null_restart;
	iter->get_next = null_get_next;
	iter->has_more = null_has_more;
	iter->set_dst = null_set_dst;

	char *memblock = iter_ptr;
	iter->dst_iter = memblock + sizeof(nb_iterator_t);
}

static uint16_t null_get_memsize(void)
{
	return 0;
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

static void null_create(void)
{
	; /* Null statement */
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
	iter->copy(iter->dst_iter, src_iter->dst_iter);
}

static void copy_handlers(nb_iterator_t *iter, const nb_iterator_t *src_iter)
{
	iter->get_memsize = src_iter->get_memsize;
	iter->init = src_iter->init;
	iter->copy = src_iter->copy;
	iter->clear = src_iter->clear;
	iter->create = src_iter->create;
	iter->clone = src_iter->clone;
	iter->destroy = src_iter->destroy;
	iter->restart = src_iter->restart;
	iter->get_next = src_iter->get_next;
	iter->has_more = src_iter->has_more;
}

inline void nb_iterator_clear(void* iter_ptr)
{
	nb_iterator_t *iter = iter_ptr;
	if (NULL != iter->dst_iter)
		iter->clear(iter->dst_iter);
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
	}
	iter->init(iter->dst_iter);
	iter->set_dst(iter->dst_iter, container->dst);
}

inline void nb_iterator_restart(nb_iterator_t *iter)
{
	iter->restart(iter->dst_iter);
}

inline const void* nb_iterator_get_next(nb_iterator_t *iter)
{
	return iter->get_next(iter->dst_iter);
}

inline bool nb_iterator_has_more(const nb_iterator_t *const iter)
{
	return iter->has_more(iter->dst_iter);
}
