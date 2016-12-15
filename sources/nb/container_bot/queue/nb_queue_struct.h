#ifndef __NB_CONTAINER_BOT_QUEUE_NB_QUEUE_STRUCT_H__
#define __NB_CONTAINER_BOT_QUEUE_NB_QUEUE_STRUCT_H__

#include "nb/memory_bot.h"
#include <stdint.h>

typedef struct {
	/* Circular linked list */
	uint32_t length;
	nb_membank_t *membank;
	nb_queue_node_t *end;
} nb_queue_t;

#endif
