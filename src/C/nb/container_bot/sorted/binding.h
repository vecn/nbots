#ifndef __NB_CONTAINER_BOT_SORTED_BINDING_H__
#define __NB_CONTAINER_BOT_SORTED_BINDING_H__

#include <stdint.h>

#include "nb/container_bot/container.h"
#include "nb/container_bot/iterator.h"

void sorted_set_handlers(nb_container_t *container);

void sorted_iterator_set_handlers(nb_iterator_t *iter);

void* sorted_do(nb_container_t *container, const char* func,
		void *data, int8_t *status);

#endif
