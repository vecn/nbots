#ifndef __NB_MEMORY_BOT_NB_MALLOC_H__
#define __NB_MEMORY_BOT_NB_MALLOC_H__

#include <stdint.h>

void *nb_malloc(uint64_t size);
void *nb_calloc(uint64_t size);
void *nb_realloc(void *mem, uint64_t new_size);
void nb_free(void *mem);

#endif
