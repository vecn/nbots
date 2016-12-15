#ifndef __NB_CONTAINER_BOT_STACK_STACK_NODE_H__
#define __NB_CONTAINER_BOT_STACK_STACK_NODE_H__

#include <stdint.h>

#include "nb/memory_bot.h"

typedef struct stack_node_s stack_node_t;

struct stack_node_s {
	stack_node_t *next;
	void *val;
};

uint32_t stack_node_get_memsize(void);
stack_node_t* stack_node_create(nb_membank_t *membank);
void stack_node_destroy(nb_membank_t *membank, stack_node_t *node);
stack_node_t* stack_node_clone(nb_membank_t *membank,
			       const stack_node_t *const node,
			       void* (*clone)(const void *const));
stack_node_t* stack_node_get_prev(const stack_node_t *const node);

#endif
