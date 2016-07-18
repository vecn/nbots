#include <stdint.h>
#include <stdlib.h>

#include "nb/memory_bot.h"

#include "stack_node.h"

uint32_t stack_node_get_memsize(void)
{
	return sizeof(stack_node_t);
}

stack_node_t* stack_node_create(nb_membank_t *membank)
{
  	return nb_membank_data_calloc(membank);
}

void stack_node_destroy(nb_membank_t *membank, stack_node_t *node)
{
	nb_membank_data_free(membank, node);
}

stack_node_t* stack_node_clone(nb_membank_t *membank,
			       const stack_node_t *const node,
			       void* (*clone)(const void *const))
{
  	stack_node_t *nc = stack_node_create(membank);
  	nc->val = clone(node->val);
  	nc->next = node->next;
  	return nc;
}

stack_node_t* stack_node_get_prev(const stack_node_t *const node)
{
	stack_node_t *iter = node->next;
	while (iter->next != node)
		iter = iter->next;
	return iter;
}
