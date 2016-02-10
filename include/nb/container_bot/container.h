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

typedef struct vcn_container_s vcn_container_t;

enum {
  VCN_CONTAINER_QUEUE,
  VCN_CONTAINER_STACK,
  VCN_CONTAINER_SORTED,
  VCN_CONTAINER_HEAP,
  VCN_CONTAINER_HASH,
  VCN_CONTAINER_NULL
};
  
vcn_container_t* vcn_container_create(int8_t id);
vcn_container_t* vcn_container_clone
			(const vcn_container_t *const container);
void vcn_container_merge(vcn_container_t *main,
			 vcn_container_t *to_clear);
void vcn_container_cast(vcn_container_t *container, int8_t new_id);
void** vcn_container_cast_to_array(vcn_container_t *container);
void vcn_container_destroy(vcn_container_t *container);
void vcn_container_clear(vcn_container_t *container);
void vcn_container_copy_to_array(const vcn_container_t *const cont_src,
				 void **array_dest);
void vcn_container_set_key_generator(vcn_container_t *container,
				     uint32_t (*key)(const void *const));
void vcn_container_set_destroyer(vcn_container_t *container,
				 void (*destroy)(void*));
void vcn_container_set_comparer(vcn_container_t *container,
				bool (*are_equal)(const void *const, 
						  const void *const));
void vcn_container_set_cloner(vcn_container_t *container,
			      void* (*clone)(const void *const));
bool vcn_container_insert(vcn_container_t *container,
			  const void *const val);
void vcn_container_insert_array(vcn_container_t *container, 
				uint32_t N, void **array);
void* vcn_container_get_first(const vcn_container_t *const container);
void* vcn_container_delete_first(vcn_container_t *container);
void* vcn_container_exist(const vcn_container_t *const container,
			  const void *const val);
void* vcn_container_delete(vcn_container_t *container,
			   const void *const val);
uint32_t vcn_container_get_length
			(const vcn_container_t *const container);
bool vcn_container_is_empty(const vcn_container_t *const container);
bool vcn_container_is_not_empty(const vcn_container_t *const container);
int8_t vcn_container_get_id(const vcn_container_t *const container);
void* vcn_container_do(vcn_container_t *container, const char* func,
		       void *data, int8_t *status);
void* vcn_container_get_dst(const vcn_container_t *const container);

#endif
