#ifndef __NB_MODEL2D_EDGE_H__
#define __NB_MODEL2D_EDGE_H__

#include <stdint.h>
#include <stdbool.h>

#include "vtx.h"

typedef struct {
	vtx_t *v1;
	vtx_t *v2;
	double length;
} edge_t;

edge_t* edge_create(void);
edge_t* edge_clone(const edge_t *const edge);
void edge_destroy(void *edge);
void edge_set_length(edge_t *edge);
double edge_get_length(const edge_t *const edge);
uint32_t edge_key_by_length(const void *const edge_ptr);
bool edge_are_equal(const void *const edge1_ptr, const void *const edge2_ptr);

#endif
