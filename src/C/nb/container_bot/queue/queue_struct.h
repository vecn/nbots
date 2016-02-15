#ifnedf __NB_CONTAINER_BOT_QUEUE_STRUCT_H__
#define __NB_CONTAINER_BOT_QUEUE_STRUCT_H__

typedef struct {
	/* Circular linked list */
	uint32_t length;
	node_t *end;
} queue_t;

#endif
