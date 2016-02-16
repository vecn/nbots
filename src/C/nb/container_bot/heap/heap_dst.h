#ifndef __NB_CONTAINER_BOT_HEAP_HEAP_DST_H__
#define __NB_CONTAINER_BOT_HEAP_HEAP_DST_H__

#include <stdbool.h>
#include <stdint.h>

uint16_t heap_get_memsize(void);
void heap_init(void *heap_ptr);
void heap_copy(void *heap_ptr, const void *src_heap_ptr,
	       void* (*clone)(const void*));
void heap_finish(void *heap_ptr,
		 void (*destroy)(void*));
void* heap_create(void);
void* heap_clone(const void *const heap_ptr,
		 void* (*clone)(const void*));
void heap_destroy(void *heap_ptr,
		  void (*destroy)(void*));
void heap_clear(void *heap_ptr,
		void (*destroy)(void*));
void heap_merge(void *heap1_ptr, void *heap2_ptr,
		uint32_t (*key)(const void*),
		int8_t (*compare)(const void*, const void*));
bool heap_insert(void *heap_ptr, const void *val,
		 uint32_t (*key)(const void*),
		 int8_t (*compare)(const void*, const void*));
void* heap_get_first(const void *const heap);
void* heap_delete_first(void *heap_ptr,
			uint32_t (*key)(const void *const));
void* heap_exist(const void *const heap_ptr, const void *val,
		 uint32_t (*key)(const void*),
		 int8_t (*compare)(const void*, const void*));
void* heap_delete(void *heap_ptr, const void *val,
		  uint32_t (*key)(const void*),
		  int8_t (*compare)(const void*, const void*));
uint32_t heap_get_length(const void *const heap_ptr);
bool heap_is_empty(const void *const heap_ptr);
bool heap_is_not_empty(const void *const heap_ptr);

#endif
