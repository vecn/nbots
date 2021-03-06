#ifndef __NB_CONTAINER_BOT_STACK_BINDING_H__
#define __NB_CONTAINER_BOT_STACK_BINDING_H__

#include <stdint.h>

#include "nb/container_bot/container.h"
#include "nb/container_bot/iterator.h"

#include "stack_dst.h"
#include "stack_iterator.h"

void stack_set_handlers(nb_container_t *container);

void stack_iterator_set_handlers(nb_iterator_t *iter);

void* stack_do(nb_container_t *container, const char* func,
	       void *data, int8_t *status);
#endif
