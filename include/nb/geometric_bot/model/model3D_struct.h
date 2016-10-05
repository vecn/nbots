#ifndef __NB_GEOMETRIC_BOT_MODEL_MODEL3D_STRUCT_H__
#define __NB_GEOMETRIC_BOT_MODEL_MODEL3D_STRUCT_H__

#include <stdint.h>

typedef struct {
	/* STereoLithography (STL) */
	uint32_t N_vtx;
	double *vtx;

	uint32_t N_face;
	double *nf;
	uint32_t *adj;
	
} nb_model3D_t;

#endif
