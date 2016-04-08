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

#include "../queue/queue_dst.h"
#include "../queue/queue_iterator.h"
#include "hash_dst.h"
#include "hash_struct.h"

static void init(hash_t *hash, float max_load_factor, uint32_t size);
static void clone_lists(hash_t *hash, const hash_t *const src_hash,
			void* (*clone)(const void *const));
static void clear_rows(hash_t *hash,
		       void (*destroy)(void*));
static void* malloc_hash(void);
static void check_realloc(hash_t *hash,
			  uint32_t (*key)(const void *const),
			  int8_t (*compare)(const void*, const void*));
static float get_load_factor(hash_t *hash);
static void realloc_rows(hash_t *hash,
			 uint32_t (*key)(const void *const),
			 int8_t (*compare)(const void*, const void*));
static void reinsert_list(hash_t *hash, void *queue_ptr,
			  uint32_t (*key)(const void *const),
			  int8_t (*compare)(const void*, const void*));
static void insert(hash_t *hash, const void *const val,
		   uint32_t (*key)(const void *const),
		   int8_t (*compare)(const void*, const void*));
static uint32_t hash_key(const void *const val, uint32_t N,
			 uint32_t (*key)(const void *const));
static void merge_rows(hash_t *merge, hash_t *hash,
		       uint32_t (*key)(const void *const),
		       int8_t (*compare)(const void*, const void*));
static void null_destroy(void *val);
static nb_container_t* get_container_from_list(const void *const list);

inline uint16_t hash_get_memsize(void)
{
	return sizeof(hash_t);
}

inline void hash_init(void *hash_ptr)
{
	init(hash_ptr, 0.85, 1000);
}

static void init(hash_t *hash, float max_load_factor, uint32_t size)
{
	hash->max_load_factor = max_load_factor;
	hash->size = size;
	hash->length = 0;
	hash->rows = calloc(size, sizeof(*(hash->rows)));
}

void hash_copy(void *hash_ptr, const void *src_hash_ptr,
	       void* (*clone)(const void *const))
{
	hash_t *hash = hash_ptr;
	const hash_t *src_hash = src_hash_ptr;
	hash->max_load_factor = src_hash->max_load_factor;
	hash->size = src_hash->size;
	hash->length = src_hash->length;
	hash->rows = calloc(hash->size, sizeof(*(hash->rows)));
	clone_lists(hash, src_hash, clone);
}

static void clone_lists(hash_t *hash, const hash_t *const src_hash,
			void* (*clone)(const void *const))
{
	for (uint32_t i = 0; i < src_hash->size; i++) {
		if (NULL != src_hash->rows[i])
			hash->rows[i] = queue_clone(src_hash->rows[i], clone);
	}
}

void hash_finish(void *hash_ptr,
		 void (*destroy)(void*))
{
	hash_clear(hash_ptr, destroy);
	hash_t* hash = hash_ptr;
	if (NULL != hash->rows) {
		free(hash->rows);
		hash->rows = NULL;
	}	
}

static void clear_rows(hash_t *hash,
		       void (*destroy)(void*))
{
	for (uint32_t i = 0; i < hash->size; i++) {
		if (NULL != hash->rows[i]) {
			queue_destroy(hash->rows[i], destroy);
			hash->rows[i] = NULL;
		}
	}
}

inline void* hash_create(void)
{
	hash_t *hash = malloc_hash();
	hash_init(hash);
	return hash;
}

static inline void* malloc_hash(void)
{
	uint16_t size = hash_get_memsize();
	return malloc(size);
}

inline void* hash_clone(const void *const hash_ptr,
			void* (*clone)(const void*))
{
	hash_t *hash = malloc_hash();
	hash_copy(hash, hash_ptr, clone);
	return hash;
}

void hash_destroy(void* hash_ptr,
		    void (*destroy)(void*))
{
	hash_t *hash = hash_ptr;
	hash_finish(hash, destroy);
	free(hash);
}

void hash_clear(void* hash_ptr,
		void (*destroy)(void*))
{
	hash_t* hash = hash_ptr;
	hash->length = 0;
	if (NULL != hash->rows)
		clear_rows(hash, destroy);
}

void hash_merge(void *hash1_ptr, void *hash2_ptr,
		uint32_t (*key)(const void*),
		int8_t (*compare)(const void*, const void*))
{
  	hash_t *hash1 = (hash_t*) hash1_ptr;
	hash_t *hash2 = (hash_t*) hash2_ptr;
	if (hash_is_not_empty(hash2)) {
		hash1->length += hash2->length;
		check_realloc(hash1, key, compare);
		merge_rows(hash1, hash2, key, compare);
		hash2->length = 0;
	}
}

static void check_realloc(hash_t *hash,
			  uint32_t (*key)(const void *const),
			  int8_t (*compare)(const void*, const void*))
{
  	float load = get_load_factor(hash);
	if (load > hash->max_load_factor)
	  realloc_rows(hash, key, compare);
}

static inline float get_load_factor(hash_t *hash)
{
	return (float) hash->length / (float) hash->size;
}

static void realloc_rows(hash_t *hash,
			 uint32_t (*key)(const void *const),
			 int8_t (*compare)(const void*, const void*))
{
	void **rows = hash->rows;
	uint32_t N = hash->size;
	hash->size = N * 2;
	hash->rows = calloc(hash->size, sizeof(*(hash->rows)));

	for (uint32_t i = 0; i < N; i++) {
		if (NULL != rows[i])
		  reinsert_list(hash, rows[i], key, compare);
	}
	free(rows);
}

static void reinsert_list(hash_t *hash, void *queue_ptr,
			  uint32_t (*key)(const void *const),
			  int8_t (*compare)(const void*, const void*))
{
	while (queue_is_not_empty(queue_ptr)) {
		void *val = queue_delete_first(queue_ptr, key);
		insert(hash, val, key, compare);
	}
	queue_destroy(queue_ptr, null_destroy);
}

static void insert(hash_t *hash, const void *const val,
		   uint32_t (*key)(const void *const),
		   int8_t (*compare)(const void*, const void*))
{
	uint32_t vkey = hash_key(val, hash->size, key);
	if (NULL == hash->rows[vkey])
		hash->rows[vkey] = queue_create();
	queue_insert(hash->rows[vkey], val, key, compare);
}

static uint32_t hash_key(const void *const val, uint32_t N,
			 uint32_t (*key)(const void *const))
{
	return key(val) % N;
}

static inline void merge_rows(hash_t *merge, hash_t *hash,
			      uint32_t (*key)(const void *const),
			      int8_t (*compare)(const void*, const void*))
{
	for (uint32_t i = 0; i < hash->size; i++) {
		if (NULL != hash->rows[i]) {
		  reinsert_list(merge, hash->rows[i], key, compare);
			hash->rows[i] = NULL;
		}
	}
}

bool hash_insert(void *hash_ptr, const void *const val,
		 uint32_t (*key)(const void *const),
		 int8_t (*compare)(const void*, const void*))
{
	hash_t* hash = hash_ptr;
	hash->length += 1;
	check_realloc(hash, key, compare);
	insert(hash, val, key, compare);
	return true;
}

void* hash_get_first(const void *const hash_ptr)
{
	const hash_t *const hash = hash_ptr;
	void *val = NULL;
	if (0 < hash->length) {
		uint32_t i = 0;
		while (NULL == hash->rows[i])
			i++;
		val = queue_get_first(hash->rows[i]);
	}
	return val;
}

void* hash_delete_first(void *hash_ptr,
			  uint32_t (*key)(const void *const))
{
	hash_t *hash = hash_ptr;
	void *val = NULL;
	if (0 < hash->length) {
		uint32_t i = 0;
		while (NULL == hash->rows[i])
			i++;

		val = queue_delete_first(hash->rows[i], key);
		hash->length -= 1;
		if (queue_is_empty(hash->rows[i])) {
			queue_destroy(hash->rows[i], null_destroy);
			hash->rows[i] = NULL;
		}
	}
	return val;
}

static inline void null_destroy(void *val)
{
	; /* Null statement */
}

void* hash_exist(const void *const hash_ptr, const void *val,
		 uint32_t (*key)(const void*),
		 int8_t (*compare)(const void*, const void*))
{
	const hash_t *const hash = hash_ptr;
	void *item = NULL;
	if (0 < hash->length) {
		uint32_t hkey = hash_key(val, hash->size, key);
		if (NULL != hash->rows[hkey])
			item = queue_exist(hash->rows[hkey], val, 
					   key, compare);
	}
	return item;
}

void* hash_delete(void *hash_ptr, const void *const val,
		  uint32_t (*key)(const void*),
		  int8_t (*compare)(const void*, const void*))
{
	hash_t *hash = hash_ptr;
	void *item = NULL;
	if (0 < hash->length) {
		uint32_t i = hash_key(val, hash->size, key);
		if (NULL != hash->rows[i]) {
			item = queue_delete(hash->rows[i], val,
					    key, compare);
			if (NULL != item)
				hash->length -= 1;
			if(queue_is_empty(hash->rows[i])){
				queue_destroy(hash->rows[i], null_destroy);
				hash->rows[i] = NULL;
			}
		}
	}
	return item;
}

inline uint32_t hash_get_length(const void *const hash_ptr)
{
	const hash_t *const hash = hash_ptr;
	return hash->length;
}

inline inline bool hash_is_empty(const void *const hash_ptr)
{
	const hash_t *const hash = hash_ptr;
	return (0 == hash->length);
}

inline inline bool hash_is_not_empty(const void *const hash_ptr)
{
	const hash_t *const hash = hash_ptr;
	return (0 < hash->length);
}

inline uint32_t hash_get_size(const void *const hash_ptr)
{
	const hash_t *const hash = hash_ptr;
	return hash->size;
}

uint32_t hash_get_N_collisions(const void *const hash_ptr)
{
	const hash_t *const hash = hash_ptr;
	uint32_t N = 0;
 	for (uint32_t i = 0; i < hash->size; i++) {
		if (NULL != hash->rows[i]) {
			if (queue_get_length(hash->rows[i]) > 1)
				N += 1;
		}
	}
	return N;
}

nb_container_t* hash_get_collisions(const void *const hash_ptr)
{
	const hash_t *const hash = hash_ptr;
	nb_container_t* cnt = nb_container_create(NB_QUEUE);
 	for (uint32_t i = 0; i < hash->size; i++) {
		void *rows = hash->rows[i];
		if (NULL != rows) {
			if (queue_get_length(rows) > 1) {
				void* collisions = 
					get_container_from_list(rows);
				nb_container_insert(cnt, collisions);
			}
		}
	}
	return cnt;
}

static nb_container_t* get_container_from_list(const void *const list)
{
	nb_container_t *items = nb_container_create(NB_QUEUE);
	void *iter = queue_iter_create();
	queue_iter_set_dst(iter, list);
	while (queue_iter_has_more(iter)) {
		const void *item = queue_iter_get_next(iter);
		nb_container_insert(items, item);
	}
	queue_iter_destroy(iter);
	return items;
}
