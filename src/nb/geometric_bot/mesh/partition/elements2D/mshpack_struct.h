#ifndef _NB_GEOMETRIC_BOT_MESH_PARTITION_ELEMENTS2D_MSHPACK_STRUCT_H_
#define _NB_GEOMETRIC_BOT_MESH_PARTITION_ELEMENTS2D_MSHPACK_STRUCT_H_

#include <stdint.h>

struct nb_mshpack_s {
	uint32_t N_elems;
	double* cen;
	double* radii;
	uint32_t* N_ngb;
	uint32_t** ngb;
};

#endif
