#include <stdlib.h>
#include "nb_queue_node.h"

nb_queue_node_t* nb_queue_node_create(void)
{
  	return calloc(1, sizeof(nb_queue_node_t));
}

void nb_queue_node_destroy(nb_queue_node_t *node, void (*destroy)(void*))
{
  	destroy(node->val);
  	free(node);
}

nb_queue_node_t *nb_queue_node_clone(const nb_queue_node_t *const node,
				     void* (*clone)(const void *const))
{
	nb_queue_node_t *nc = nb_queue_node_create();
  	nc->val = clone(node->val);
  	nc->next = node->next;
  	return nc;
}

nb_queue_node_t *nb_queue_node_get_prev(const nb_queue_node_t *const node)
{
	nb_queue_node_t *iter = node->next;
	while (iter->next != node)
		iter = iter->next;
	return iter;
}
