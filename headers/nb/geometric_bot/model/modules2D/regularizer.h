#ifndef __NB_GEOMETRIC_BOT_MODEL_MODULES2D_REGULARIZER_H__
#define __NB_GEOMETRIC_BOT_MODEL_MODULES2D_REGULARIZER_H__

#include <stdint.h>
#include "nb/geometric_bot/model/model2D.h"
	
int nb_model_regularize(nb_model_t* model, double lambda,
			 uint32_t N_fixed_vertices,
			 uint32_t* fixed_vertices);

#endif
