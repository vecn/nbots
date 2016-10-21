#ifndef _NB_GEOMETRIC_BOT_MESH_MESH2D_ELEMENTS2D_MSH3TRG_STRUCT_H_
#define _NB_GEOMETRIC_BOT_MESH_MESH2D_ELEMENTS2D_MSH3TRG_STRUCT_H_

#include <stdint.h>

struct nb_msh3trg_s {
	uint32_t N_nod;
	double *nod;

	uint32_t N_edg;
	uint32_t *edg;

	uint32_t N_elems;
	uint32_t *adj;
	uint32_t *ngb;

	uint32_t N_vtx;
	uint32_t *vtx;

	uint32_t N_sgm;
	uint32_t *N_nod_x_sgm;
	uint32_t **nod_x_sgm;
};

struct nb_msh3trg_t {
	uint32_t N_nod;
	double *nod;

	uint32_t N_edg;
	uint32_t *edg;

	uint32_t N_elems;
	uint32_t *adj;
	uint32_t *ngb;

	uint32_t N_vtx;
	uint32_t *vtx;

	uint32_t N_sgm;
	uint32_t *N_nod_x_sgm;
	uint32_t **nod_x_sgm;
};

#endif
