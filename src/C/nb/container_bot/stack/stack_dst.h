#ifndef __NB_CONTAINER_BOT_STACK_DST_H__
#define __NB_CONTAINER_BOT_STACK_DST_H__

#include <stdbool.h>
#include <stdint.h>

uint16_t stack_get_memsize(void);
void stack_init(void *stack_ptr);
void stack_copy(void *stack_ptr, const void *src_stack_ptr,
		void* (*clone)(const void*));
void stack_clear(void *stack_ptr,
		 void (*destroy)(void*));
void* stack_create(void);
void* stack_clone(const void *const stack_ptr,
		 void* (*clone)(const void*));
void stack_destroy(void *stack_ptr,
		  void (*destroy)(void*));
void stack_merge(void *list1_ptr, void *list2_ptr,
		 uint32_t (*key)(const void*),
		 int8_t (*compare)(const void*, const void*));
void stack_clear(void *stack_ptr,
		void (*destroy)(void*));
bool stack_insert(void *stack_ptr, const void *val,
		  uint32_t (*key)(const void*),
		  int8_t (*compare)(const void*, const void*));
void* stack_get_first(const void *const list);
void* stack_delete_first(void *stack_ptr,
			uint32_t (*key)(const void*));
void* stack_exist(const void *const stack_ptr, const void *const val,
		 uint32_t (*key)(const void*),
		 int8_t (*compare)(const void*, const void*));
void* stack_delete(void *stack_ptr, const void *val,
		  uint32_t (*key)(const void*),
		  int8_t (*compare)(const void*, const void*));
uint32_t stack_get_length(const void *const stack_ptr);
bool stack_is_empty(const void *const stack_ptr);
bool stack_is_not_empty(const void *const stack_ptr);

#endif
