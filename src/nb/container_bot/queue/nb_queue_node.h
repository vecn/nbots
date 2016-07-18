#ifndef __NB_CONTAINER_BOT_QUEUE_NB_QUEUE_NODE_H__
#define __NB_CONTAINER_BOT_QUEUE_NB_QUEUE_NODE_H__

#include <stdint.h>

#include "nb/memory_bot.h"

typedef struct nb_queue_node_s nb_queue_node_t;

struct nb_queue_node_s {
	nb_queue_node_t *next;
	void *val;
};

uint32_t nb_queue_node_get_memsize(void);
nb_queue_node_t* nb_queue_node_create(nb_membank_t *membank);
void nb_queue_node_destroy(nb_membank_t *membank,
			   nb_queue_node_t *node);
nb_queue_node_t* nb_queue_node_clone(nb_membank_t *membank,
				     const nb_queue_node_t *const node,
				     void* (*clone)(const void *const));
nb_queue_node_t* node_get_prev(const nb_queue_node_t *const node);

#endif
