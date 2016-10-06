#ifndef __NB_MEMORY_BOT_NB_ALLOCATE_MEM_H__
#define __NB_MEMORY_BOT_NB_ALLOCATE_MEM_H__

#include <stdint.h>
#include <alloca.h>

#define MAX_STACK_MEM_ON_SOFT_MALLOCS 1024 /* 1 kB */

#define NB_SOFT_MALLOC(size)				\
	((MAX_STACK_MEM_ON_SOFT_MALLOCS < (size)) ?	\
	 nb_allocate_mem((size)):alloca((size)))
#define NB_SOFT_FREE(size,ptr)				\
	((MAX_STACK_MEM_ON_SOFT_MALLOCS < (size)) ?	\
	 nb_free_mem((ptr)):((ptr) = 0))

void *nb_allocate_mem(uint64_t size);
void *nb_allocate_zero_mem(uint64_t size);
void *nb_reallocte_mem(void *mem, uint64_t new_size);
void nb_free_mem(void *mem);

#endif
