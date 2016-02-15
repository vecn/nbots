#ifndef __NB_CONTAINER_BOT_QUEUE_QUEUE_STRUCT_H__
#define __NB_CONTAINER_BOT_QUEUE_QUEUE_STRUCT_H__

#include <stdint.h>

typedef struct {
	/* Circular linked list */
	uint32_t length;
	node_t *end;
} queue_t;

#endif
