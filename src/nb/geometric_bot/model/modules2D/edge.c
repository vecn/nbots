#include <stdlib.h>

#include "nb/geometric_bot/utils2D.h"
#include "edge.h"

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
	bool eq_e1v1_e2v1 = edge1->v1 == edge2->v1;
	bool eq_e1v2_e2v2 = edge1->v2 == edge2->v2;
	bool eq_e1v1_e2v2 = edge1->v1 == edge2->v2;
	bool eq_e1v2_e2v1 = edge1->v2 == edge2->v1;
	int8_t out;
	if ((eq_e1v1_e2v1 && eq_e1v2_e2v2) ||
	    (eq_e1v1_e2v2 && eq_e1v2_e2v1)) {
		out = 0;
	} else {
		if (edge1_ptr < edge2_ptr)
			out = -1;
		else
			out = 1;
	}
	return out;
}
