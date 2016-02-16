#ifndef __NB_CONTAINER_BOT_SORTED_AVL_DST_H__
#define __NB_CONTAINER_BOT_SORTED_AVL_DST_H__

#include <stdbool.h>
#include <stdint.h>

uint16_t avl_get_memsize(void);
void avl_init(void *avl_ptr);
void avl_copy(void *avl_ptr, const void *src_avl_ptr,
	      void* (*clone)(const void*));
void avl_clear(void *avl_ptr,
	       void (*destroy)(void*));
void* avl_create(void);
void* avl_clone(const void *const avl_ptr,
		void* (*clone)(const void*));
void avl_destroy(void *avl_ptr,
		 void (*destroy)(void*));
void avl_merge(void *avl1_ptr, void *avl2_ptr,
	       uint32_t (*key)(const void*),
	       int8_t (*compare)(const void*, const void*));
void avl_clear(void *avl_ptr,
	       void (*destroy)(void*));
bool avl_insert(void *avl_ptr, const void *val,
		uint32_t (*key)(const void*),
		int8_t (*compare)(const void*, const void*));
void* avl_get_first(const void *const avl_ptr);
void* avl_delete_first(void *avl_ptr,
		       uint32_t (*key)(const void*));
void* avl_delete_last(void *avl_ptr,
		      uint32_t (*key)(const void*));
void* avl_exist(const void *const avl_ptr, const void *val,
		uint32_t (*key)(const void*),
		int8_t (*compare)(const void*, const void*));
void* avl_delete(void *avl_ptr, const void *const val,
		 uint32_t (*key)(const void *const),
		 int8_t (*compare)(const void*, const void*));
uint32_t avl_get_length(const void *const avl_ptr);
bool avl_is_empty(const void *const avl_ptr);
bool avl_is_not_empty(const void *const avl_ptr);

#endif
