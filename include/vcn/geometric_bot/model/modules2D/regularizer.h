#ifndef __VCN_MODEL2D_REGULARIZER_H__
#define __VCN_MODEL2D_REGULARIZER_H__

#include <stdint.h>
#include "vcn/geometric_bot/model/model2D.h"

#ifdef __cplusplus
extern "C" {
#endif
	
	int vcn_model_regularize(vcn_model_t* model, double lambda,
				 uint32_t N_fixed_vertices,
				 uint32_t* fixed_vertices);

#ifdef __cplusplus
}
#endif

#endif
