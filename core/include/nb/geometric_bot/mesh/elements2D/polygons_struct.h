#ifndef __NB_MESH2D_POLYGONS_STRUCT_H__
#define __NB_MESH2D_POLYGONS_STRUCT_H__

#include <stdint.h>

typedef struct {
	uint32_t N_nod;
	double *nod;      /* Nodes coordinates concatenated */

	uint32_t N_edg;
	uint32_t *edg;

	uint32_t N_elems;
	double *cen;      /* Centroids */
	uint16_t *N_adj;
	uint32_t **adj;   /* Connectivity matrix */
	uint32_t **ngb;   /* Neighbours */

	uint32_t N_vtx;
	uint32_t *elem_vtx; /* Ids of elems corresponding to input vtx */

	uint32_t N_sgm;
	/* Number of nodes forming the input segment */
	uint32_t *N_nod_x_sgm;
	/* Sequence of nodal ids forming the input segments */
	uint32_t **nod_x_sgm;
} nb_mshpoly_t;

#endif
