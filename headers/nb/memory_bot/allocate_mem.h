#ifndef __NB_MEMORY_BOT_NB_ALLOCATE_MEM_H__
#define __NB_MEMORY_BOT_NB_ALLOCATE_MEM_H__

#include <stdint.h>

#include "nb/sysinfo.h"

#if defined(CC_Microsoft) || (defined(CC_GNU) && defined(OS_Windows))
	#include <malloc.h>
        #define nb_allocate_on_stack(size) _alloca(size)
#else

	#include <alloca.h>
	#define nb_allocate_on_stack(size) alloca(size)
#endif

#define NB_MAX_STACK_MEM_ON_SOFT_ALLOCATION 1024 /* 1 kB */

#define nb_soft_allocate_mem(size)				\
	((NB_MAX_STACK_MEM_ON_SOFT_ALLOCATION < (size)) ?       \
	 nb_allocate_mem(size):nb_allocate_on_stack(size))

#define nb_soft_free_mem(size, ptr)				\
	((NB_MAX_STACK_MEM_ON_SOFT_ALLOCATION < (size)) ?	\
	 nb_free_mem(ptr):((ptr) = 0))

void *nb_allocate_mem(uint64_t size);
void *nb_allocate_zero_mem(uint64_t size);
void *nb_reallocate_mem(void *mem, uint64_t new_size);
void nb_free_mem(void *mem);

#endif
