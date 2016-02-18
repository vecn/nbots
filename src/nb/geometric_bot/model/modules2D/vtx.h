#ifndef __NB_MODEL2D_VTX_H__
#define __NB_MODEL2D_VTX_H__

#include <stdint.h>
#include <stdbool.h>

typedef struct {
	uint32_t id;
	double x[2];
  	bool initial_input;
} vtx_t;

vtx_t* vtx_create(void);
void vtx_destroy(void *vtx);
void vtx_set_id(vtx_t *vtx, uint32_t id);
uint32_t vtx_get_id(const vtx_t *const vtx);
bool vtx_is_not_initial(const vtx_t *const vtx);
void vtx_set_as_initial(vtx_t *vtx);
uint32_t vtx_hash_key(const void *const  vertex);
int8_t vtx_compare(const void *const vtxA,
		   const void *const vtxB);

#endif
