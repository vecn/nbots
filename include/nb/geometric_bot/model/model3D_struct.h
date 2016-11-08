#ifndef __NB_GEOMETRIC_BOT_MODEL_MODEL3D_STRUCT_H__
#define __NB_GEOMETRIC_BOT_MODEL_MODEL3D_STRUCT_H__

#include <stdint.h>

typedef struct {
	/* STereoLithography (STL) */
	uint32_t N_vtx;                   /*Number of vertex*/
	double *vtx;                      /*Train of vertexes*/

	uint32_t N_face;                  /*Number of faces*/
	double *nf;                       /*Trainf of normal for face*/
	uint32_t *adj;                    /*Train of node index by face*/
	
} nb_model3D_t;


#endif
