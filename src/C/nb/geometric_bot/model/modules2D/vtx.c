#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>

#include "nb/geometric_bot/utils2D.h"
#include "vtx.h"

inline vtx_t* vtx_create(void)
{
	return calloc(1, sizeof(vtx_t));
}

inline void vtx_destroy(void *vtx)
{
	free(vtx);
}

inline void vtx_set_id(vtx_t *vtx, uint32_t id)
{
	vtx->id = id;
}

inline uint32_t vtx_get_id(const vtx_t *const vtx)
{
	return vtx->id;
}

inline bool vtx_is_not_initial(const vtx_t *const vtx)
{
	return !vtx->initial_input;
}

inline void vtx_set_as_initial(vtx_t *vtx)
{
	vtx->initial_input = true;
}

inline uint32_t vtx_hash_key(const void *const  vertex)
{
	const vtx_t *const vtx = vertex;
	return (uint32_t)((int)(vtx->x[0] * 73856093) ^
			  (int)(vtx->x[1] * 19349663));
}

int8_t vtx_compare(const void *const vtxA, const void *const vtxB)
{
	const vtx_t *const vA = vtxA;
	const vtx_t *const vB = vtxB;
	return (vcn_utils2D_get_dist2(vA->x, vB->x) < NB_GEOMETRIC_TOL) ?
	  0:1;
}
