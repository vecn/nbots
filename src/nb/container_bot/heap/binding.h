#ifndef __NB_CONTAINER_BOT_HEAP_BINDING_H__
#define __NB_CONTAINER_BOT_HEAP_BINDING_H__

#include <stdint.h>

#include "nb/container_bot/container.h"
#include "nb/container_bot/iterator.h"

#include "heap_dst.h"
#include "heap_iterator.h"

void heap_set_handlers(nb_container_t *container);

void heap_iterator_set_handlers(nb_iterator_t *iter);

void* heap_do(nb_container_t *container, const char* func,
		void *data, int8_t *status);

#endif
