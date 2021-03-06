#ifndef __NB_CONTAINER_BOT_QUEUE_NB_QUEUE_BINDING_H__
#define __NB_CONTAINER_BOT_QUEUE_NB_QUEUE_BINDING_H__

#include <stdint.h>

#include "nb/container_bot/container.h"
#include "nb/container_bot/iterator.h"

#include "nb_queue_dst.h"
#include "nb_queue_iterator.h"

void nb_queue_set_handlers(nb_container_t *container);

void nb_queue_iterator_set_handlers(nb_iterator_t *iter);

void* nb_queue_do(nb_container_t *container, const char* func,
		  void *data, int8_t *status);

#endif
