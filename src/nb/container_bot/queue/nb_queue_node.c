#include <stdlib.h>
#include "nb_queue_node.h"

uint32_t nb_queue_node_get_memsize(void)
{
	return sizeof(nb_queue_node_t);
}

nb_queue_node_t* nb_queue_node_create(nb_membank_t *membank)
{
  	return nb_membank_calloc(membank);
}

void nb_queue_node_destroy(nb_membank_t *membank,
			   nb_queue_node_t *node)
{
	nb_membank_free(membank, node);
}

nb_queue_node_t *nb_queue_node_clone(nb_membank_t *membank,
				     const nb_queue_node_t *const node,
				     void* (*clone)(const void *const))
{
	nb_queue_node_t *nc = nb_queue_node_create(membank);
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
