#ifndef __NB_MEMORY_BOT_NB_ALLOCATE_MEM_H__
#define __NB_MEMORY_BOT_NB_ALLOCATE_MEM_H__

#include <stdint.h>

#include "nb/sysinfo.h"

#ifdef __FRAMAC__
        #define nb_allocate_on_stack(size) malloc(size)
        #define nb_free_on_stack(ptr) free(ptr)
#else
        #if defined(CC_Microsoft) || (defined(CC_GNU) && defined(OS_Windows))
	        #include <malloc.h>
                #define nb_allocate_on_stack(size) _alloca(size)
                #define nb_free_on_stack(ptr) ((ptr)=0)
        #else
	        #include <alloca.h>
	        #define nb_allocate_on_stack(size) alloca(size)
                #define nb_free_on_stack(ptr) ((ptr)=0)
        #endif
#endif

#define NB_MAX_STACK_MEM_ON_SOFT_ALLOCATION 1024 /* 1 kB */

#define nb_soft_allocate_mem(size)				\
	((NB_MAX_STACK_MEM_ON_SOFT_ALLOCATION < (size)) ?       \
	 nb_allocate_mem(size):nb_allocate_on_stack(size))

#define nb_soft_free_mem(size, ptr)				\
	((NB_MAX_STACK_MEM_ON_SOFT_ALLOCATION < (size)) ?	\
	 nb_free_mem(ptr):nb_free_on_stack(ptr))

void *nb_allocate_mem(uint64_t size);
void *nb_allocate_zero_mem(uint64_t size);
void *nb_reallocate_mem(void *mem, uint64_t new_size);
void nb_free_mem(void *mem);

#endif
