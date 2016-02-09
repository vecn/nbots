#ifndef __VCN_HTABLE_DST_H__
#define __VCN_HTABLE_DST_H__

#include <stdbool.h>
#include <stdint.h>

#include "vcn/container_bot/container.h"

void* htable_create(void);
void* htable_clone(const void *const htable_ptr,
		   void* (*clone)(const void *const));
void htable_merge(void *htable1_ptr, void *htable2_ptr,
		  uint32_t (*key)(const void *const));
void htable_destroy(void *htable_ptr,
		    void (*destroy)(void*));
void htable_clear(void *htable_ptr,
		  void (*destroy)(void*));
bool htable_insert(void *htable_ptr, const void *const val,
		   uint32_t (*key)(const void *const));
void* htable_get_first(const void *const htable_ptr);
void* htable_delete_first(void *htable_ptr,
			  uint32_t (*key)(const void *const));
void* htable_exist(const void *const htable_ptr, const void *const val,
		   uint32_t (*key)(const void *const),
		   bool (*are_equal)(const void *const, const void *const));
void* htable_delete(void *htable_ptr, const void *const val,
		    uint32_t (*key)(const void *const),
		    bool (*are_equal)(const void *const, const void *const));
uint32_t htable_get_length(const void *const htable_ptr);
bool htable_is_empty(const void *const htable_ptr);
bool htable_is_not_empty(const void *const htable_ptr);
const void* htable_get_iterator_start(const void *const htable_ptr);
uint32_t htable_get_size(const void *const htable_ptr);
uint32_t htable_get_N_collisions(const void *const htable_ptr);
vcn_container_t* htable_get_collisions(const void *const htable_ptr);

#endif
