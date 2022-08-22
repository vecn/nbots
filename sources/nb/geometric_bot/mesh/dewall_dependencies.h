#ifndef __NB_GEOMETRIC_BOT_MESH_DEWALL_DEPENDENCIES_H__
#define __NB_GEOMETRIC_BOT_MESH_DEWALL_DEPENDENCIES_H__

#include <stdint.h>
#include <stdbool.h>

typedef struct {
	void* (*allocate)(uint64_t);
	void* (*allocate_zero)(uint64_t);
	void (*free)(void*);
} interface_heap_t;

typedef void* afl_t; /* Active Face List (Set) */

typedef struct {
	uint32_t (*size)(void);
	void (*init)(afl_t*, uint32_t (*keygen)(const void *const));
	void (*finish)(afl_t*);
	bool (*is_empty)(const afl_t *const);
	void (*insert)(afl_t*, const void *const);
	void* (*delete_any)(afl_t*);	
	void* (*delete)(afl_t*, const void *const);
} interface_afl_t;

typedef struct {
	interface_heap_t mem;
	interface_afl_t afl;
} interface_t;

interface_t* module(void);

#endif
