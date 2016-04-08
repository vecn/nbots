#ifndef __NB_CONTAINER_BOT_QUEUE_DST_H__
#define __NB_CONTAINER_BOT_QUEUE_DST_H__

#include <stdbool.h>
#include <stdint.h>

uint16_t queue_get_memsize(void);
void queue_init(void *queue_ptr);
void queue_copy(void *queue_ptr, const void *src_queue_ptr,
		void* (*clone)(const void*));
void queue_finish(void *queue_ptr,
		  void (*destroy)(void*));
void* queue_create(void);
void* queue_clone(const void *const queue_ptr,
		 void* (*clone)(const void*));
void queue_destroy(void *queue_ptr,
		  void (*destroy)(void*));
void queue_clear(void *queue_ptr,
		 void (*destroy)(void*));
void queue_merge(void *list1_ptr, void *list2_ptr,
		 uint32_t (*key)(const void*),
	       int8_t (*compare)(const void*, const void*));
void queue_clear(void *queue_ptr,
		void (*destroy)(void*));
bool queue_insert_first(void *queue_ptr, const void *val,
			uint32_t (*key)(const void*),
			int8_t (*compare)(const void*, const void*));
bool queue_insert(void *queue_ptr, const void *val,
		  uint32_t (*key)(const void*),
		  int8_t (*compare)(const void*, const void*));
void* queue_get_first(const void *const list);
void* queue_delete_first(void *queue_ptr,
			uint32_t (*key)(const void*));
void* queue_exist(const void *const queue_ptr, const void *const val,
		 uint32_t (*key)(const void*),
		 int8_t (*compare)(const void*, const void*));
void* queue_delete(void *queue_ptr, const void *val,
		  uint32_t (*key)(const void*),
		  int8_t (*compare)(const void*, const void*));
uint32_t queue_get_length(const void *const queue_ptr);
bool queue_is_empty(const void *const queue_ptr);
bool queue_is_not_empty(const void *const queue_ptr);

#endif
