#ifndef __NB_HEAP_DST_H__
#define __NB_HEAP_DST_H__

#include <stdbool.h>
#include <stdint.h>

void* heap_create(void);
void* heap_clone(const void *const heap_ptr,
		 void* (*clone)(const void *const));
void heap_merge(void *heap1_ptr, void *heap2_ptr,
		uint32_t (*key)(const void *const));
void heap_destroy(void *heap_ptr,
		  void (*destroy)(void*));
void heap_clear(void *heap_ptr,
		void (*destroy)(void*));
bool heap_insert(void *heap_ptr, const void *const val,
		 uint32_t (*key)(const void *const));
void* heap_get_first(const void *const heap);
void* heap_delete_first(void *heap_ptr,
			uint32_t (*key)(const void *const));
void* heap_exist(const void *const heap_ptr, const void *const val,
		 uint32_t (*key)(const void *const),
		 bool (*are_equal)(const void *const, const void *const));
void* heap_delete(void *heap_ptr, const void *const val,
		  uint32_t (*key)(const void *const),
		  bool (*are_equal)(const void *const, const void *const));
uint32_t heap_get_length(const void *const heap_ptr);
bool heap_is_empty(const void *const heap_ptr);
bool heap_is_not_empty(const void *const heap_ptr);
const void* heap_get_iterator_start(const void *const heap_ptr);

#endif
