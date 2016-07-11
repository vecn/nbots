#include <stdlib.h>
#include <stdint.h>

#include "nb/memory_bot/nb_malloc.h"

void *nb_malloc(uint64_t size)
{
	return malloc(size);
}

void *nb_calloc(uint64_t size)
{
	return calloc(size, 1);
}

void *nb_realloc(void *mem, uint64_t new_size)
{
	return realloc(mem, new_size);
}

void nb_free(void *mem)
{
	return free(mem);
}
