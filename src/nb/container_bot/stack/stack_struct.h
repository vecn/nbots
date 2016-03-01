#ifndef __NB_CONTAINER_BOT_STACK_STACK_STRUCT_H__
#define __NB_CONTAINER_BOT_STACK_STACK_STRUCT_H__

#include <stdint.h>

#include "stack_node.h"

typedef struct {
	/* Circular linked list */
	uint32_t length;
	stack_node_t *end;
} stack_t;

#endif
