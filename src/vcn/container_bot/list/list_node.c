#include <stdlib.h>
#include "list_node.h"

inline node_t* node_create(void)
{
  	return calloc(1, sizeof(node_t));
}

inline void node_destroy(node_t *node, void (*destroy)(void*))
{
  	destroy(node->val);
  	free(node);
}

node_t* node_clone(const node_t *const node,
		   void* (*clone)(const void *const))
{
  	node_t *nc = node_create();
  	nc->val = clone(node->val);
  	nc->next = node->next;
  	return nc;
}

node_t* node_get_prev(const node_t *const node)
{
	node_t *iter = node->next;
	while (iter->next != node)
		iter = iter->next;
	return iter;
}
