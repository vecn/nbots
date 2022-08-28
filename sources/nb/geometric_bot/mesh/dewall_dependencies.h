#ifndef __NB_GEOMETRIC_BOT_MESH_DEWALL_DEPENDENCIES_H__
#define __NB_GEOMETRIC_BOT_MESH_DEWALL_DEPENDENCIES_H__

#include <stdint.h>
#include <stdbool.h>

typedef struct {
	void* (*allocate)(uint64_t);
	void* (*allocate_zero)(uint64_t);
	void (*free)(void*);
} interface_heap_t;

typedef void afl_t; /* Active Face List (Set) */

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
	void (*init)(afl_iterator_t*, const afl_t *const);
	void (*finish)(afl_iterator_t*);
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

typedef void mesh_t; /* Output */
typedef msh_trg_t trg_t;// TODO: typedef void trg_t;
typedef msh_vtx_t vtx_t;// TODO: typedef void vtx_t;
typedef struct {
	bool (*is_v3_intersecting_any_edge)(const mesh_t *const,
					    const vtx_t *const v1,
					    const vtx_t *const v2,
					    const vtx_t *const v3);
	trg_t* (*new_triangle)(mesh_t*, /* Counter-clockwise order */
			       const vtx_t *const v1,
			       const vtx_t *const v2,
			       const vtx_t *const v3);
} interface_mesh_t;

typedef struct {
	interface_heap_t mem; // TODO: Remove and pass membank instead
	interface_afl_t afl;
	interface_afl_iterator_t afl_iter;
	interface_queue_t queue;
	interface_mesh_t mesh;  /* Output handler */
} interface_t;

/* WARNING: 
 *   module() must not allocate memory,
 *   it should return a reference to static structure.
 */
interface_t* module(void);

#endif
