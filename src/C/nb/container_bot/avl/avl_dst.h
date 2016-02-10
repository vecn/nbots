#ifndef __NB_AVL_DST_H__
#define __NB_AVL_DST_H__

#include <stdbool.h>
#include <stdint.h>

void* avl_create(void);
void* avl_clone(const void *const avl_ptr,
		void* (*clone)(const void *const));
void avl_merge(void *avl1_ptr, void *avl2_ptr,
	       uint32_t (*key)(const void *const));
void avl_destroy(void *avl_ptr,
		 void (*destroy)(void*));
void avl_clear(void *avl_ptr,
	       void (*destroy)(void*));
bool avl_insert(void *avl_ptr, const void *const val,
		uint32_t (*key)(const void *const));
void* avl_get_first(const void *const avl_ptr);
void* avl_delete_first(void *avl_ptr,
		       uint32_t (*key)(const void *const));
void* avl_delete_last(void *avl_ptr,
		      uint32_t (*key)(const void *const));
void* avl_exist(const void *const avl_ptr, const void *const val,
		uint32_t (*key)(const void *const),
		bool (*are_equal)(const void *const, const void *const));
void* avl_delete(void *avl_ptr, const void *const val,
		 uint32_t (*key)(const void *const),
		 bool (*are_equal)(const void *const, const void *const));
uint32_t avl_get_length(const void *const avl_ptr);
bool avl_is_empty(const void *const avl_ptr);
bool avl_is_not_empty(const void *const avl_ptr);
const void* avl_get_iterator_start(const void *const avl_ptr);

#endif
