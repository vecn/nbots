#ifndef __NB_CONTAINER_BOT_HASH_BINDING_H__
#define __NB_CONTAINER_BOT_HASH_BINDING_H__

#include <stdint.h>

#include "nb/container_bot/container.h"
#include "nb/container_bot/iterator.h"

void hash_set_handlers(nb_container_t *container);

void hash_iterator_set_handlers(nb_iterator_t *iter);

void* hash_do(nb_container_t *container, const char* func,
		void *data, int8_t *status);

#endif
