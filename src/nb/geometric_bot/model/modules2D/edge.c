#include <stdlib.h>

#include "nb/geometric_bot/utils2D.h"
#include "edge.h"

#define ABS(a) (((a) > 0)?(a):(-(a)))
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

inline edge_t* edge_create(void)
{
	return calloc(1, sizeof(edge_t));
}

edge_t* edge_clone(const edge_t *const edge)
{
	edge_t *ec = edge_create();
	ec->length = edge->length;
	ec->v1 = edge->v1;
	ec->v2 = edge->v2;
	return ec;
}

inline void edge_destroy(void *edge)
{
	free(edge);
}

inline void edge_set_length(edge_t *edge)
{
	edge->length = vcn_utils2D_get_dist(edge->v1->x, edge->v2->x);
}

inline double edge_get_length(const edge_t *const edge)
{
	return edge->length;
}

inline uint32_t edge_key_by_length(const void *const edge_ptr)
{
	return (uint32_t)(1e8 * edge_get_length(edge_ptr));
}

int8_t edge_compare(const void *const edge1_ptr, const void *const edge2_ptr)
{
	const edge_t *const edge1 = edge1_ptr;
	const edge_t *const edge2 = edge2_ptr;
	int8_t out;

	bool eq_e1v1_e2v1 = edge1->v1 == edge2->v1;
	bool eq_e1v2_e2v2 = edge1->v2 == edge2->v2;
	bool eq_e1v1_e2v2 = edge1->v1 == edge2->v2;
	bool eq_e1v2_e2v1 = edge1->v2 == edge2->v1;
	if ((eq_e1v1_e2v1 && eq_e1v2_e2v2) ||
	    (eq_e1v1_e2v2 && eq_e1v2_e2v1)) {
		out = 0;
		goto EXIT;
	}

	int32_t xv1e1 = (int32_t)(1000 * edge1->v1->x[0]);
	int32_t xv2e1 = (int32_t)(1000 * edge1->v2->x[0]);
	int32_t abs1 = ABS(xv1e1 - xv2e1);
	int32_t xv1e2 = (int32_t)(1000 * edge2->v1->x[0]);
	int32_t xv2e2 = (int32_t)(1000 * edge2->v2->x[0]);
	int32_t abs2 = ABS(xv1e2 - xv2e2);
	if (abs1 < abs2) {
		out = 1;
		goto EXIT;
	} else if (abs1 > abs2) {
		out = -1;
		goto EXIT;
	}
	
	int32_t yv1e1 = (int32_t)(1000 * edge1->v1->x[1]);
	int32_t yv2e1 = (int32_t)(1000 * edge1->v2->x[1]);
	abs1 = ABS(yv1e1 - yv2e1);
	int32_t yv1e2 = (int32_t)(1000 * edge2->v1->x[1]);
	int32_t yv2e2 = (int32_t)(1000 * edge2->v2->x[1]);
	abs2 = ABS(yv1e2 - yv2e2);
	if (abs1 < abs2) {
		out = 1;
		goto EXIT;
	} else if (abs1 > abs2) {
		out = -1;
		goto EXIT;
	}

	int32_t c1 = MIN(xv1e1, xv2e1);
	int32_t c2 = MIN(xv1e2, xv2e2);
	if (c1 < c2) {
		out = 1;
		goto EXIT;
	} else if (c1 > c2) {
		out = -1;
		goto EXIT;
	}
	
	c1 = MAX(xv1e1, xv2e1);
	c2 = MAX(xv1e2, xv2e2);
	if (c1 < c2) {
		out = 1;
		goto EXIT;
	} else if (c1 > c2) {
		out = -1;
		goto EXIT;
	}

	c1 = MIN(yv1e1, yv2e1);
	c2 = MIN(yv1e2, yv2e2);
	if (c1 < c2) {
		out = 1;
		goto EXIT;
	} else if (c1 > c2) {
		out = -1;
		goto EXIT;
	}
	
	c1 = MAX(yv1e1, yv2e1);
	c2 = MAX(yv1e2, yv2e2);
	if (c1 < c2) {
		out = 1;
		goto EXIT;
	} else if (c1 > c2) {
		out = -1;
		goto EXIT;
	}
	out = 0;
EXIT:
	return out;
}
