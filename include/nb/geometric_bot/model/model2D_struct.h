#ifndef __NB_GEOMETRIC_BOT_MODEL_MODEL2D_STRUCT_H__
#define __NB_GEOMETRIC_BOT_MODEL_MODEL2D_STRUCT_H__

#include <stdint.h>

typedef struct {
	/* Planar Straight Line Graph (PSLG) */
	uint32_t N;              /* Number of vertices */
	uint32_t M;              /* Number of edges */
	uint32_t H;              /* Number of holes */
	double *vertex;          /* Vertex coordinates sorted by id */
	uint32_t *edge;          /* Pairs of vertex ids forming an edge 
				  * sorted by id */
	double *holes;       /* Coordinates inside holes */
} nb_model_t;

void nb_model_alloc_vertices(nb_model_t *model);
void nb_model_alloc_edges(nb_model_t *model);
void nb_model_alloc_holes(nb_model_t *model);

#endif
