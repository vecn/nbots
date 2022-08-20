#include "nb/memory_bot.h"

#include "dewall_dependencies.h"

static interface_t implementation = {
	nb_allocate_mem,
	nb_allocate_zero_mem,
	nb_free_mem
};

interface_t* module(void)
{
	return &implementation;
}
