#ifndef __NB_MEMORY_BOT_NB_MALLOC_H__
#define __NB_MEMORY_BOT_NB_MALLOC_H__

#include <stdint.h>
#include <alloca.h>

#define MAX_STACK_MEM_ON_SOFT_MALLOCS 1024 /* 1 kB */

#define NB_SOFT_MALLOC(size)				\
	((MAX_STACK_MEM_ON_SOFT_MALLOCS < (size)) ?	\
	 nb_malloc((size)):alloca((size)))
#define NB_SOFT_FREE(size,ptr)				\
	((MAX_STACK_MEM_ON_SOFT_MALLOCS < (size)) ?	\
	 nb_free((ptr)):((ptr) = 0))

void *nb_malloc(uint64_t size);
void *nb_calloc(uint64_t size);
void *nb_realloc(void *mem, uint64_t new_size);
void nb_free(void *mem);

#endif
