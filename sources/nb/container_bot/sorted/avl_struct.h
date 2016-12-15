#ifndef __NB_CONTAINER_BOT_SORTED_AVL_STRUCT_H__
#define __NB_CONTAINER_BOT_SORTED_AVL_STRUCT_H__

#include "nb/memory_bot.h"

typedef struct {
	uint32_t length;
	nb_membank_t *membank;
	tree_t* root;
} avl_t;

#endif
