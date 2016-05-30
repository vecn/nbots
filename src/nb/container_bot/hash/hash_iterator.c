/******************************************************************************
 *   Hash-Table iterator                                                      *
 *   2011-2015 Victor Eduardo Cardoso Nungaray                                *
 *   Twitter: @victore_cardoso                                                *
 *   email: victorc@cimat.mx                                                  *
 ******************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "../queue/queue_dst.h"
#include "../queue/queue_iterator.h"
#include "hash_struct.h"
#include "hash_iterator.h"

typedef struct {
	/* Iterator for the hash table */
	void ** ht_rows;  /* Pointer to the rows of the hash table */
	uint32_t ht_size; /* Size of rows */
	uint32_t i;             /* Index of the hash table */
	void* row_iter;         /* Row iterator */
} iter_t;

static void set_next_row_iterator(iter_t *iter);
static void* malloc_hash_iter(void);

inline uint16_t hash_iter_get_memsize(void)
{
	return sizeof(iter_t);
}

inline void hash_iter_init(void *iter_ptr)
{
	memset(iter_ptr, 0, sizeof(iter_t));
}

void hash_iter_copy(void *iter_ptr, const void *src_iter_ptr)
{
	iter_t *iter = iter_ptr;
	const iter_t *src_iter = src_iter_ptr;
	iter->ht_rows = src_iter->ht_rows;
	iter->ht_size = src_iter->ht_size;
	iter->i = src_iter->i;
	iter->row_iter = queue_iter_clone(src_iter->row_iter);
}

void hash_iter_finish(void *iter_ptr)
{
	hash_iter_clear(iter_ptr);
}

void* hash_iter_create(void)
{
	void *iter = malloc_hash_iter();
	hash_iter_init(iter);
	return iter;
}

static inline void* malloc_hash_iter(void)
{
	uint16_t size = hash_iter_get_memsize();
	return malloc(size);
}

void* hash_iter_clone(const void *iter_ptr)
{
	void *iter = malloc_hash_iter();
	hash_iter_copy(iter, iter_ptr);
	return iter;
}

void hash_iter_destroy(void *iter_ptr)
{
	hash_iter_finish(iter_ptr);
	free(iter_ptr);
}

void hash_iter_clear(void *iter_ptr)
{
	iter_t *iter = iter_ptr;
	if (NULL != iter->row_iter)
		queue_iter_destroy(iter->row_iter);
	memset(iter, 0, sizeof(iter_t));
}

void hash_iter_set_dst(void *iter_ptr, const void *hash_ptr)
{
	iter_t *iter = iter_ptr;
	if (NULL != hash_ptr) {
		const hash_t* hash = hash_ptr;
		iter->ht_rows = hash->rows; 
		iter->ht_size = hash->size;
		if (NULL == iter->row_iter)
			iter->row_iter = queue_iter_create();
		hash_iter_restart(iter);
	}
}

inline void hash_iter_restart(void* iter_ptr)
{
	iter_t *iter = iter_ptr;
	iter->i = 0;
	if (NULL != iter->row_iter)	
		set_next_row_iterator(iter);
}

static void set_next_row_iterator(iter_t *iter)
{
	while (iter->i < iter->ht_size) {
		if (NULL != iter->ht_rows[iter->i]) {
			queue_iter_set_dst(iter->row_iter,
					  iter->ht_rows[iter->i]);
			break;
		}
		iter->i += 1;
	}
	if (iter->i >= iter->ht_size)
		queue_iter_set_dst(iter->row_iter, NULL);

}

const void* hash_iter_get_next(void *iter_ptr)
{
	iter_t *iter = iter_ptr;
	const void *val = NULL;
	if (NULL != iter->row_iter)
		val = queue_iter_get_next(iter->row_iter);
	return val;
}

bool hash_iter_has_more(const void *const iter_ptr)
{
	iter_t *iter = (iter_t*) iter_ptr;
	bool has_more = false;
	if (NULL != iter->row_iter) {
		has_more = queue_iter_has_more(iter->row_iter);
		if (!has_more) {
			iter->i += 1;
			set_next_row_iterator(iter);
			has_more = queue_iter_has_more(iter->row_iter);
		}
	}
	return has_more;
}
