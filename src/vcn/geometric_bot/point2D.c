#include <stdlib.h>

#include "vcn/geometric_bot/utils2D.h"
#include "vcn/geometric_bot/point2D.h"

inline vcn_point2D_t* vcn_point2D_create(void)
{
	return calloc(1, sizeof(vcn_point2D_t));
}

inline void vcn_point2D_destroy(vcn_point2D_t* point)
{
	free(point);
}

inline bool vcn_point2D_are_equal(const void *const p1_ptr,
				  const void *const p2_ptr)
{
	const vcn_point2D_t *const p1 = p1_ptr;
	const vcn_point2D_t *const p2 = p2_ptr;
	return vcn_utils2D_get_dist(p1->x, p2->x) < VCN_GEOMETRIC_TOL_POW2;
}
