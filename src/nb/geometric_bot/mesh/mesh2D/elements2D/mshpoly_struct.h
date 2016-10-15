#ifndef _NB_GEOMETRIC_BOT_MESH_MESH2D_ELEMENTS2D_MSHPOLY_STRUCT_H_
#define _NB_GEOMETRIC_BOT_MESH_MESH2D_ELEMENTS2D_MSHPOLY_STRUCT_H_

#include <stdint.h>

struct nb_mshpoly_s {
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
	uint32_t *vtx; /* Ids of elems corresponding to input vtx.
			     * 'N_elems' for input vtx on the boundary. */

	uint32_t N_sgm;
	/* Number of nodes forming the input segment */
	uint32_t *N_nod_x_sgm;
	/* Sequence of nodal ids forming the input segments */
	uint32_t **nod_x_sgm;
};

#endif
