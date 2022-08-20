#ifndef __NB_GEOMETRIC_BOT_MESH_DEWALL_DEPENDENCIES_H__
#define __NB_GEOMETRIC_BOT_MESH_DEWALL_DEPENDENCIES_H__

#include <stdint.h>

typedef struct {
	void* (*allocate)(uint64_t);
	void* (*allocate_zero)(uint64_t);
	void (*free)(void*);
} interface_t;

interface_t* module(void);

#endif
