#ifndef __NB_CONTAINER_BOT_STACK_STACK_NODE_H__
#define __NB_CONTAINER_BOT_STACK_STACK_NODE_H__

typedef struct stack_node_s stack_node_t;

struct stack_node_s {
	stack_node_t *next;
	void *val;
};

stack_node_t* stack_node_create(void);
void stack_node_destroy(stack_node_t *node);
stack_node_t* stack_node_clone(const stack_node_t *const node,
		    void* (*clone)(const void *const));
stack_node_t* stack_node_get_prev(const stack_node_t *const node);

#endif
