#include <stdlib.h>
#include <stdint.h>
#include <alloca.h>

#include "nb/memory_bot/allocate_mem.h"

void *nb_allocate_mem(uint64_t size)
{
	return malloc(size);
}

void *nb_allocate_zero_mem(uint64_t size)
{
	return calloc(size, 1);
}

void *nb_reallocate_mem(void *mem, uint64_t new_size)
{
	return realloc(mem, new_size);
}

void nb_free_mem(void *mem)
{
	return free(mem);
}
