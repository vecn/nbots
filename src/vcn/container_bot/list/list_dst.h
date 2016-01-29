#ifndef __VCN_LIST_DST_H__
#define __VCN_LIST_DST_H__

#include <stdbool.h>
#include <stdint.h>

void* list_create(void);
void* list_clone(const void *const list_ptr,
		 void* (*clone)(const void *const));
void list_merge(void *list1_ptr, void *list2_ptr,
		 uint32_t (*key)(const void *const));
void list_destroy(void *list_ptr,
		  void (*destroy)(void*));
void list_clear(void *list_ptr,
		void (*destroy)(void*));
bool list_insert_first(void *list_ptr, const void *const val,
		       uint32_t (*key)(const void *const));
bool list_insert_last(void *list_ptr, const void *const val,
		      uint32_t (*key)(const void *const));
void* list_get_first(const void *const list);
void* list_delete_first(void *list_ptr,
			uint32_t (*key)(const void *const));
void* list_exist(const void *const list_ptr, const void *const val,
		 uint32_t (*key)(const void *const),
		 bool (*are_equal)(const void *const, const void *const));
void* list_delete(void *list_ptr, const void *const val,
		  uint32_t (*key)(const void *const),
		  bool (*are_equal)(const void *const, const void *const));
uint32_t list_get_length(const void *const list_ptr);
bool list_is_empty(const void *const list_ptr);
bool list_is_not_empty(const void *const list_ptr);
const void* list_get_iterator_start(const void *const list_ptr);

#endif
