/******************************************************************************
 *   Hash Table DST: Hash Table                                               *
 *   2011-2015 Victor Eduardo Cardoso Nungaray                                *
 *   Twitter: @victore_cardoso                                                *
 *   email: victorc@cimat.mx                                                  *
 ******************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>

#include "nb/container_bot/container.h"

#include "../list/list_dst.h"
#include "../list/list_iterator.h"
#include "htable_dst.h"

typedef struct {
	/* Hash Table */
	float max_load_factor;
	uint32_t length;
	uint32_t size;
	void** rows;/* Pointers to lists */
} htable_t;

static htable_t* create(float max_load_factor, uint32_t size);
static void clone_lists(htable_t *cloned, const htable_t *const htable,
			void* (*clone)(const void *const));
static void check_realloc(htable_t *htable,
			  uint32_t (*key)(const void *const));
static float get_load_factor(htable_t *htable);
static void realloc_rows(htable_t *htable,
			 uint32_t (*key)(const void *const));
static void reinsert_list(htable_t *htable, void *list_ptr,
			  uint32_t (*key)(const void *const));
static void insert(htable_t *htable, const void *const val,
		   uint32_t (*key)(const void *const));
static uint32_t hash_key(const void *const val, uint32_t N,
			 uint32_t (*key)(const void *const));
static void merge_rows(htable_t *merge, htable_t *htable,
		       uint32_t (*key)(const void *const));
static void clear_rows(htable_t *htable,
		       void (*destroy)(void*));
static void null_destroy(void *val);
static vcn_container_t* get_container_from_list(const void *const list);

inline void* htable_create(void)
{
  	return create(0.85, 1000);
}

static htable_t* create(float max_load_factor, uint32_t size)
{
  	htable_t *htable = calloc(1, sizeof(*htable));
	htable->max_load_factor = max_load_factor;
	htable->size = size;
	htable->rows = calloc(size, sizeof(*(htable->rows)));
	return htable;
}

void* htable_clone(const void *const htable_ptr,
		   void* (*clone)(const void *const))
{
  	const htable_t *const htable = htable_ptr;
	htable_t *hc = create(htable->max_load_factor, htable->size);
	hc->length = htable->length;
	clone_lists(hc, htable, clone);
	return hc;
}

static void clone_lists(htable_t *cloned, const htable_t *const htable,
			void* (*clone)(const void *const))
{
	for (uint32_t i = 0; i < htable->size; i++) {
		if (NULL != htable->rows[i])
			cloned->rows[i] = list_clone(htable->rows[i], clone);
	}
}

void htable_merge(void *htable1_ptr, void *htable2_ptr,
		  uint32_t (*key)(const void *const))
{
  	htable_t *htable1 = (htable_t*) htable1_ptr;
	htable_t *htable2 = (htable_t*) htable2_ptr;
	if (htable_is_not_empty(htable2)) {
		htable1->length += htable2->length;
		check_realloc(htable1, key);
		merge_rows(htable1, htable2, key);
		htable2->length = 0;
	}
}

static void check_realloc(htable_t *htable,
			  uint32_t (*key)(const void *const))
{
  	float load = get_load_factor(htable);
	if (load > htable->max_load_factor)
		realloc_rows(htable, key);
}

static inline float get_load_factor(htable_t *htable)
{
	return (float) htable->length / (float) htable->size;
}

static void realloc_rows(htable_t *htable,
			 uint32_t (*key)(const void *const))
{
	void **rows = htable->rows;
	uint32_t N = htable->size;
	htable->size = N * 2;
	htable->rows = calloc(htable->size, sizeof(*(htable->rows)));

	for (uint32_t i = 0; i < N; i++) {
		if (NULL != rows[i])
			reinsert_list(htable, rows[i], key);
	}
	free(rows);
}

static void reinsert_list(htable_t *htable, void *list_ptr,
			  uint32_t (*key)(const void *const))
{
	while (list_is_not_empty(list_ptr)) {
		void *val = list_delete_first(list_ptr, key);
		insert(htable, val, key);
	}
	list_destroy(list_ptr, null_destroy);
}

static void insert(htable_t *htable, const void *const val,
		   uint32_t (*key)(const void *const))
{
	uint32_t vkey = hash_key(val, htable->size, key);
	if (NULL == htable->rows[vkey])
		htable->rows[vkey] = list_create();
	list_insert_last(htable->rows[vkey], val, key);
}

static uint32_t hash_key(const void *const val, uint32_t N,
			 uint32_t (*key)(const void *const))
{
	return key(val) % N;
}

static inline void merge_rows(htable_t *merge, htable_t *htable,
			      uint32_t (*key)(const void *const))
{
	for (uint32_t i = 0; i < htable->size; i++) {
		if (NULL != htable->rows[i]) {
			reinsert_list(merge, htable->rows[i], key);
			htable->rows[i] = NULL;
		}
	}
}

void htable_destroy(void* htable_ptr,
		    void (*destroy)(void*))
{
	htable_t *htable = htable_ptr;
	htable_clear(htable, destroy);
	free(htable->rows);
	free(htable);
}

void htable_clear(void* htable_ptr,
		  void (*destroy)(void*))
{
	htable_t* htable = htable_ptr;
	clear_rows(htable, destroy);
	htable->length = 0;  
}

static void clear_rows(htable_t *htable,
		       void (*destroy)(void*))
{
	for (uint32_t i = 0; i < htable->size; i++) {
		if (NULL != htable->rows[i]) {
			list_destroy(htable->rows[i], destroy);
			htable->rows[i] = NULL;
		}
	}
}

bool htable_insert(void *htable_ptr, const void *const val,
		   uint32_t (*key)(const void *const))
{
	htable_t* htable = htable_ptr;
	htable->length += 1;
	check_realloc(htable, key);
	insert(htable, val, key);
	return true;
}

void* htable_get_first(const void *const htable_ptr)
{
	const htable_t *const htable = htable_ptr;
	void *val = NULL;
	if (0 < htable->length) {
		uint32_t i = 0;
		while (NULL == htable->rows[i])
			i++;
		val = list_get_first(htable->rows[i]);
	}
	return val;
}

void* htable_delete_first(void *htable_ptr,
			  uint32_t (*key)(const void *const))
{
	htable_t *htable = htable_ptr;
	void *val = NULL;
	if (0 < htable->length) {
		uint32_t i = 0;
		while (NULL == htable->rows[i])
			i++;

		val = list_delete_first(htable->rows[i], key);
		htable->length -= 1;
		if (list_is_empty(htable->rows[i])) {
			list_destroy(htable->rows[i], null_destroy);
			htable->rows[i] = NULL;
		}
	}
	return val;
}

static inline void null_destroy(void *val)
{
	; /* Null statement */
}

void* htable_exist(const void *const htable_ptr, const void *const val,
		   uint32_t (*key)(const void *const),
		   bool (*are_equal)(const void *const, const void *const))
{
	const htable_t *const htable = htable_ptr;
	void *item = NULL;
	if (0 < htable->length) {
		uint32_t hkey = hash_key(val, htable->size, key);
		if (NULL != htable->rows[hkey])
			item = list_exist(htable->rows[hkey], val, 
					  key, are_equal);
	}
	return item;
}

void* htable_delete(void *htable_ptr, const void *const val,
		    uint32_t (*key)(const void *const),
		    bool (*are_equal)(const void *const, const void *const))
{
	htable_t *htable = htable_ptr;
	void *item = NULL;
	if (0 < htable->length) {
		uint32_t i = hash_key(val, htable->size, key);
		if (NULL != htable->rows[i]) {
			item = list_delete(htable->rows[i], val,
					   key, are_equal);
			if (NULL != item)
				htable->length -= 1;
			if(list_is_empty(htable->rows[i])){
				list_destroy(htable->rows[i], null_destroy);
				htable->rows[i] = NULL;
			}
		}
	}
	return item;
}

inline uint32_t htable_get_length(const void *const htable_ptr)
{
	const htable_t *const htable = htable_ptr;
	return htable->length;
}

inline inline bool htable_is_empty(const void *const htable_ptr)
{
	const htable_t *const htable = htable_ptr;
	return (0 == htable->length);
}

inline inline bool htable_is_not_empty(const void *const htable_ptr)
{
	const htable_t *const htable = htable_ptr;
	return (0 < htable->length);
}

inline const void* htable_get_iterator_start(const void *const htable_ptr)
{
	const htable_t *const htable = htable_ptr;
	return htable->rows;
}

inline uint32_t htable_get_size(const void *const htable_ptr)
{
	const htable_t *const htable = htable_ptr;
	return htable->size;
}

uint32_t htable_get_N_collisions(const void *const htable_ptr)
{
	const htable_t *const htable = htable_ptr;
	uint32_t N = 0;
 	for (uint32_t i = 0; i < htable->size; i++) {
		if (NULL != htable->rows[i]) {
			if (list_get_length(htable->rows[i]) > 1)
				N += 1;
		}
	}
	return N;
}

vcn_container_t* htable_get_collisions(const void *const htable_ptr)
{
	const htable_t *const htable = htable_ptr;
	vcn_container_t* cnt = vcn_container_create(VCN_CONTAINER_QUEUE);
 	for (uint32_t i = 0; i < htable->size; i++) {
		void *rows = htable->rows[i];
		if (NULL != rows) {
			if (list_get_length(rows) > 1) {
				void* collisions = 
					get_container_from_list(rows);
				vcn_container_insert(cnt, collisions);
			}
		}
	}
	return cnt;
}

static vcn_container_t* get_container_from_list(const void *const list)
{
	vcn_container_t *items = vcn_container_create(VCN_CONTAINER_QUEUE);
	void *iter = list_iter_create();
	list_iter_set_dst(iter, list);
	while (list_iter_has_more(iter)) {
		const void *item = list_iter_get_next(iter);
		vcn_container_insert(items, item);
	}
	list_iter_destroy(iter);
	return items;
}
