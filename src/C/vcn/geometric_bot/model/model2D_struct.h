#ifndef __VCN_GEOMETRIC_BOT_MODEL_MODEL2D_STRUCT_H__
#define __VCN_GEOMETRIC_BOT_MODEL_MODEL2D_STRUCT_H__

#include <stdint.h>
#include "vcn/geometric_bot/model/model2D.h"

struct vcn_model_s {
	/* Planar Straight Line Graph (PSLG) */
	uint32_t N;              /* Number of vertices */
	uint32_t M;              /* Number of edges */
	uint32_t H;              /* Number of holes */
	double *vertex;      /* Vertex coordinates sorted by id */
	uint32_t *edge;          /* Pairs of vertex ids forming an edge 
			      * sorted by id */
	double *holes;       /* Coordinates inside holes */
};

void model_alloc_vertices(vcn_model_t *model);
void model_alloc_edges(vcn_model_t *model);
void model_alloc_holes(vcn_model_t *model);

#endif
