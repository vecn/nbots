#ifndef __NB_GEOMETRIC_BOT_MESH_ELEMENTS2D_QUAD_STRUCT_H__
#define __NB_GEOMETRIC_BOT_MESH_ELEMENTS2D_QUAD_STRUCT_H__

#include <stdint.h>

typedef struct {
	uint32_t N_nod;
	double *nod;      /* Nodes coordinates concatenated */

	uint32_t N_edg;
	uint32_t *edg;

	uint32_t N_elems;
	int8_t *type;     /* Quad if 0, Trg otherwise. */
	uint32_t *adj;    /* Connectivity matrix (4 x N_elems) */
	uint32_t *ngb;    /* Quad-neighbours (4 x N_elems) */

	uint32_t N_vtx;
	uint32_t *vtx; /* Ids of vtx corresponding to input vtx */

	uint32_t N_sgm;
	/* Number of nodes forming the input segment */
	uint32_t *N_nod_x_sgm;
	/* Sequence of nodal ids forming the input segments */
	uint32_t** nod_x_sgm;
} nb_mshquad_t;


#endif
