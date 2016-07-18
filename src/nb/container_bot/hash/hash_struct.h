#ifndef __CONTAINER_BOT_HASH_HASH_STRUCT_H__
#define __CONTAINER_BOT_HASH_HASH_STRUCT_H__

#include "nb/memory_bot.h"

#include <stdint.h>

typedef struct {
	/* Hash Table */
	float max_load_factor;
	uint32_t length;
	uint32_t size;
	nb_membank_t *membank;
	void** rows;/* Pointers to lists */
} hash_t;

#endif
