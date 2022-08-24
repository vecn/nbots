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


typedef void* afl_iterator_t;  /* AFL Iterator */

typedef struct {
	uint32_t (*size)(void);
	void (*init)(afl_iterator_t*);
	void (*finish)(afl_iterator_t*);
	void (*set_afl)(afl_iterator_t*, const afl_t *const);
	const void* (*get_next)(afl_iterator_t*);
	bool (*has_more)(const afl_iterator_t *const);
} interface_afl_iterator_t;

typedef void *queue_t;
typedef struct {
	uint32_t (*size)(void);
	void (*init)(queue_t*);
	void (*finish)(queue_t*);
	void (*add)(queue_t*, const void *const);
	void* (*poll)(queue_t*);
	bool (*is_empty)(const queue_t *const);

} interface_queue_t;

typedef struct {
	interface_heap_t mem;
	interface_afl_t afl;
	interface_afl_iterator_t afl_iter;
	interface_queue_t queue;
} interface_t;

/* WARNING: 
 *   module() must not allocate memory,
 *   it should return a reference to static structure.
 */
interface_t* module(void);

#endif
