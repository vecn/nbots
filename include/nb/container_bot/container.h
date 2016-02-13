/**
 * @file conainer.h
 * @brief  Generic container.
 * @author Victor Eduardo Cardoso Nungaray
 * @n victorc@@cimat.mx
 * @n @@victore_cardoso
 */

#ifndef __NB_CONTAINER_BOT_CONTAINER_H__
#define __NB_CONTAINER_BOT_CONTAINER_H__

#include <stdbool.h>
#include <stdint.h>

typedef struct nb_container_s nb_container_t;

typedef enum {
  NB_QUEUE,
  NB_STACK,
  NB_SORTED,
  NB_HEAP,
  NB_HASH,
  NB_NULL
} nb_container_type;

uint16_t nb_container_get_memsize(nb_container_type type);
void nb_container_init(void *container_ptr, nb_container_type type);
void nb_container_copy(void *container_ptr, const void *src_container_ptr);
void nb_container_clear(void* container_ptr);
void* nb_container_create(nb_container_type type);
void* nb_container_clone(const void *container_ptr);
void nb_container_destroy(void *container_ptr);

void nb_container_merge(nb_container_t *main,
			 nb_container_t *to_clear);
void nb_container_cast(nb_container_t *container, int8_t new_id);
void** nb_container_cast_to_array(nb_container_t *container);
void nb_container_clear(nb_container_t *container);
void nb_container_copy_to_array(const nb_container_t *const cont_src,
				 void **array_dest);
void nb_container_set_key_generator(nb_container_t *container,
				     uint32_t (*key)(const void *const));
void nb_container_set_destroyer(nb_container_t *container,
				 void (*destroy)(void*));
void nb_container_set_comparer(nb_container_t *container,
				bool (*are_equal)(const void *const, 
						  const void *const));
void nb_container_set_cloner(nb_container_t *container,
			      void* (*clone)(const void *const));
bool nb_container_insert(nb_container_t *container,
			  const void *const val);
void nb_container_insert_array(nb_container_t *container, 
				uint32_t N, void **array);
void* nb_container_get_first(const nb_container_t *const container);
void* nb_container_delete_first(nb_container_t *container);
void* nb_container_exist(const nb_container_t *const container,
			  const void *const val);
void* nb_container_delete(nb_container_t *container,
			   const void *const val);
uint32_t nb_container_get_length
			(const nb_container_t *const container);
bool nb_container_is_empty(const nb_container_t *const container);
bool nb_container_is_not_empty(const nb_container_t *const container);
int8_t nb_container_get_id(const nb_container_t *const container);
void* nb_container_do(nb_container_t *container, const char* func,
		       void *data, int8_t *status);
void* nb_container_get_dst(const nb_container_t *const container);

#endif
