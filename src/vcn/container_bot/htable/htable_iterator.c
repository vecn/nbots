/******************************************************************************
 *   Hash-Table iterator                                                      *
 *   2011-2015 Victor Eduardo Cardoso Nungaray                                *
 *   Twitter: @victore_cardoso                                                *
 *   email: victorc@cimat.mx                                                  *
 ******************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "../list/list_dst.h"
#include "../list/list_iterator.h"
#include "htable_dst.h"
#include "htable_iterator.h"

typedef struct {
	/* Iterator for the hash table */
	void ** ht_rows;  /* Pointer to the rows of the hash table */
	uint32_t ht_size; /* Size of rows */
	uint32_t i;             /* Index of the hash table */
	void* row_iter;         /* Row iterator */
} iter_t;

static void set_next_row_iterator(iter_t *iter);

inline void* htable_iter_create(void)
{
	return calloc(1, sizeof(iter_t));
}

void htable_iter_set_dst(void *iter_ptr, const void *const htable_ptr)
{
	iter_t *iter = iter_ptr;
	if (NULL != htable_ptr) {
		iter->ht_rows = (void**) htable_get_iterator_start(htable_ptr); 
		iter->ht_size = htable_get_size(htable_ptr);
		if (NULL == iter->row_iter)
			iter->row_iter = list_iter_create();
		htable_iter_restart(iter);
	}
}

void* htable_iter_clone(const void *const iter_ptr)
{
	const iter_t *const restrict iter = iter_ptr;
	iter_t* clone = calloc(1, sizeof(*clone));
	clone->ht_rows = iter->ht_rows;
	clone->ht_size = iter->ht_size;
	clone->i = iter->i;
	clone->row_iter = list_iter_clone(iter->row_iter);
	return clone;
}

inline void htable_iter_destroy(void *iter_ptr)
{
	iter_t *iter = iter_ptr;
	if (NULL != iter->row_iter)
		list_iter_destroy(iter->row_iter);
	free(iter);
}

void htable_iter_restart(void* iter_ptr)
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
			list_iter_set_dst(iter->row_iter,
					  iter->ht_rows[iter->i]);
			break;
		}
		iter->i += 1;
	}
	if (iter->i >= iter->ht_size)
		list_iter_set_dst(iter->row_iter, NULL);

}

const void* htable_iter_get_next(void *iter_ptr)
{
	iter_t *iter = iter_ptr;
	const void *val = NULL;
	if (NULL != iter->row_iter)
		val = list_iter_get_next(iter->row_iter);
	return val;
}

bool htable_iter_has_more(const void *const iter_ptr)
{
	iter_t *iter = (iter_t*) iter_ptr;
	bool has_more = false;
	if (NULL != iter->row_iter) {
		has_more = list_iter_has_more(iter->row_iter);
		if (!has_more) {
			iter->i += 1;
			set_next_row_iterator(iter);
			has_more = list_iter_has_more(iter->row_iter);
		}
	}
	return has_more;
}
