#ifndef __NB_GEOMETRIC_BOT_MESH_DEWALL_DEPENDENCIES_H__
#define __NB_GEOMETRIC_BOT_MESH_DEWALL_DEPENDENCIES_H__

#include <stdint.h>
#include <stdbool.h>

typedef struct {
	void* (*allocate)(uint64_t);
	void* (*allocate_zero)(uint64_t);
	void (*free)(void*);
} interface_heap_t;

typedef void* hash_table_t;

typedef struct {
	uint32_t (*size)(void);
	void (*init)(hash_table_t*);
	void (*finish)(hash_table_t*);
	void (*set_keygen)(hash_table_t*, uint32_t (*k)(const void *const));
	bool (*is_empty)(const hash_table_t *const);
	void (*insert)(hash_table_t*, const void *const);
	void* (*delete_any)(hash_table_t*);	
	void* (*delete)(hash_table_t*, const void *const);
} interface_hash_table_t;

typedef struct {
	interface_heap_t mem;
	interface_hash_table_t hash_table;
} interface_t;

interface_t* module(void);

#endif
