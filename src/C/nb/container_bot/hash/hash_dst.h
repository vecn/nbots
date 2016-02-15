#ifndef __NB_CONTAINER_BOT_HASH_HASH_DST_H__
#define __NB_CONTAINER_BOT_HASH_HASH_DST_H__

#include <stdbool.h>
#include <stdint.h>

#include "nb/container_bot/container.h"

uint16_t hash_get_memsize(void);
void hash_init(void *hash_ptr);
void hash_copy(void *hash_ptr, const void *src_hash_ptr,
	       void* (*clone)(const void *const));
void hash_clear(void *hash_ptr,
		void (*destroy)(void*));
void* hash_create(void);
void* hash_clone(const void *const hash_ptr,
		 void* (*clone)(const void *const));
void hash_destroy(void *hash_ptr,
		  void (*destroy)(void*));
void hash_merge(void *hash1_ptr, void *hash2_ptr,
		uint32_t (*key)(const void *const));
bool hash_insert(void *hash_ptr, const void *const val,
		 uint32_t (*key)(const void *const));
void* hash_get_first(const void *const hash_ptr);
void* hash_delete_first(void *hash_ptr,
			uint32_t (*key)(const void *const));
void* hash_exist(const void *const hash_ptr, const void *val,
		 uint32_t (*key)(const void*),
		 int8_t (*compare)(const void*, const void*));
void* hash_delete(void *hash_ptr, const void *val,
		  uint32_t (*key)(const void*),
		  int8_t (*compare)(const void*, const void*));
uint32_t hash_get_length(const void *const hash_ptr);
bool hash_is_empty(const void *const hash_ptr);
bool hash_is_not_empty(const void *const hash_ptr);
uint32_t hash_get_size(const void *const hash_ptr);
uint32_t hash_get_N_collisions(const void *const hash_ptr);
nb_container_t* hash_get_collisions(const void *const hash_ptr);

#endif
