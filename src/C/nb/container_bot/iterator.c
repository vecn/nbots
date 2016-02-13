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
	bool (*has_more)(const void *const);
};

static void* null_clone(const void *const iter);
static void null_destroy(void *iter);
static void null_restart(void *iter);
static const void* null_get_next(void *iter);
static bool null_has_more(const void *const iter);

static void iterator_set_queue(nb_iterator_t *iter, const void *const queue);
static void iterator_set_stack(nb_iterator_t *iter, const void *const stack);
static void iterator_set_avl(nb_iterator_t *iter, const void *const avl);
static void iterator_set_htable(nb_iterator_t *iter, const void *const htable);
static void iterator_set_heap(nb_iterator_t *iter, const void *const heap);

void nb_iterator_init(nb_iterator_t *iter)
{
	iter->clone = null_clone;
	iter->destroy = null_destroy;
	iter->restart = null_restart;
	iter->get_next = null_get_next;
	iter->has_more = null_has_more;
}

nb_iterator_t* nb_iterator_create(void)
{
	nb_iterator_t *iter = calloc(1, sizeof(*iter));
	nb_iterator_init(iter);
	return iter;
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

void nb_iterator_set_container(nb_iterator_t *iter,
				const nb_container_t *const container)
{
	void *dst = nb_container_get_dst(container);
	switch (nb_container_get_id(container)) {
	case NB_QUEUE:
		iterator_set_queue(iter, dst);
		break;
	case NB_STACK:
		iterator_set_stack(iter, dst);
		break;
	case NB_SORTED:
		iterator_set_avl(iter, dst);
		break;
	case NB_HEAP:
		iterator_set_heap(iter, dst);
		break;
	case NB_HASH:
		iterator_set_htable(iter, dst);
		break;
	default:
		; /* Null statement */
	}
}

static void iterator_set_queue(nb_iterator_t *iter, const void *const queue)
{
	if (NULL != iter->dst_iter)
		iter->destroy(iter->dst_iter);

	iter->dst_iter = list_iter_create();
	list_iter_set_dst(iter->dst_iter, queue);

	iter->clone = list_iter_clone;
	iter->destroy = list_iter_destroy;
	iter->restart = list_iter_restart;
	iter->get_next = list_iter_get_next;
	iter->has_more = list_iter_has_more;
}

static void iterator_set_stack(nb_iterator_t *iter, const void *const stack)
{
	if (NULL != iter->dst_iter)
		iter->destroy(iter->dst_iter);

	iter->dst_iter = list_iter_create();
	list_iter_set_dst(iter->dst_iter, stack);

	iter->clone = list_iter_clone;
	iter->destroy = list_iter_destroy;
	iter->restart = list_iter_restart;
	iter->get_next = list_iter_get_next;
	iter->has_more = list_iter_has_more;
}

static void iterator_set_avl(nb_iterator_t *iter, const void *const avl)
{
	if (NULL != iter->dst_iter)
		iter->destroy(iter->dst_iter);

	iter->dst_iter = avl_iter_create();
	avl_iter_set_dst(iter->dst_iter, avl);

	iter->clone = avl_iter_clone;
	iter->destroy = avl_iter_destroy;
	iter->restart = avl_iter_restart;
	iter->get_next = avl_iter_get_next;
	iter->has_more = avl_iter_has_more;
}

static void iterator_set_htable(nb_iterator_t *iter, const void *const htable)
{
	if (NULL != iter->dst_iter)
		iter->destroy(iter->dst_iter);

	iter->dst_iter = htable_iter_create();
	htable_iter_set_dst(iter->dst_iter, htable);

	iter->clone = htable_iter_clone;
	iter->destroy = htable_iter_destroy;
	iter->restart = htable_iter_restart;
	iter->get_next = htable_iter_get_next;
	iter->has_more = htable_iter_has_more;
}

static void iterator_set_heap(nb_iterator_t *iter, const void *const heap)
{
	if (NULL != iter->dst_iter)
		iter->destroy(iter->dst_iter);

	iter->dst_iter = heap_iter_create();
	heap_iter_set_dst(iter->dst_iter, heap);

	iter->clone = heap_iter_clone;
	iter->destroy = heap_iter_destroy;
	iter->restart = heap_iter_restart;
	iter->get_next = heap_iter_get_next;
	iter->has_more = heap_iter_has_more;
}

nb_iterator_t* nb_iterator_clone(const nb_iterator_t *const restrict iter)
{
	nb_iterator_t *clone = malloc(sizeof(*clone));
	clone->dst_iter = iter->clone(iter->dst_iter);

	clone->clone = iter->clone;
	clone->destroy = iter->destroy;
	clone->restart = iter->restart;
	clone->get_next = iter->get_next;
	clone->has_more = iter->has_more;
	return clone;
}

inline void nb_iterator_destroy(nb_iterator_t *iter)
{
	if (NULL != iter->dst_iter)
		iter->destroy(iter->dst_iter);
	free(iter);
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
