#ifndef __NB_CONTAINER_BOT_QUEUE_NB_QUEUE_NODE_H__
#define __NB_CONTAINER_BOT_QUEUE_NB_QUEUE_NODE_H__

typedef struct nb_queue_node_s nb_queue_node_t;

struct nb_queue_node_s {
	nb_queue_node_t *next;
	void *val;
};

nb_queue_node_t* nb_queue_node_create(void);
void nb_queue_node_destroy(nb_queue_node_t *node, void (*destroy)(void*));
nb_queue_node_t* nb_queue_node_clone(const nb_queue_node_t *const node,
				     void* (*clone)(const void *const));
nb_queue_node_t* node_get_prev(const nb_queue_node_t *const node);

#endif
