#ifndef __NB_CONTAINER_BOT_STACK_STACK_STRUCT_H__
#define __NB_CONTAINER_BOT_STACK_STACK_STRUCT_H__

#include <stdint.h>

#include "nb/memory_bot.h"

#include "stack_node.h"

typedef struct {
	/* Circular linked list */
	uint32_t length;
	nb_membank_t *membank;
	stack_node_t *end;
} nb_stack_t;

#endif
