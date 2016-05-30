#ifndef __NB_CONTAINER_BOT_QUEUE_QUEUE_NODE_H__
#define __NB_CONTAINER_BOT_QUEUE_QUEUE_NODE_H__

typedef struct node_s node_t;

struct node_s {
	node_t *next;
	void *val;
};

node_t* node_create(void);
void node_destroy(node_t *node, void (*destroy)(void*));
node_t* node_clone(const node_t *const node,
		    void* (*clone)(const void *const));
node_t* node_get_prev(const node_t *const node);

#endif
