#ifndef __NB_CONTAINER_BOT_HEAP_HEAP_STRUCT_H__
#define __NB_CONTAINER_BOT_HEAP_HEAP_STRUCT_H__

#include <stdint.h>

#include "nb/memory_bot.h"

typedef struct {
	/* Pairing heap (Half-ordered binary tree) */
	uint32_t length;
	nb_membank_t *membank;
	htree_t* root;
} heap_t;

#endif
