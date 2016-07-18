#include <stdlib.h>
#include "stack_node.h"

inline stack_node_t* stack_node_create(void)
{
  	return calloc(1, sizeof(stack_node_t));
}

inline void stack_node_destroy(stack_node_t *node)
{
  	free(node);
}

stack_node_t* stack_node_clone(const stack_node_t *const node,
		   void* (*clone)(const void *const))
{
  	stack_node_t *nc = stack_node_create();
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
