#ifndef __VCN_MODEL2D_BLENDER_H__
#define __VCN_MODEL2D_BLENDER_H__

#include <stdbool.h>
#include <stdint.h>

#include "vcn/geometric_bot/mesh/elements2D/triangles_struct.h"
#include "vcn/geometric_bot/model/model2D.h"
#include "vcn/graph_bot.h"

#ifdef __cplusplus
extern "C" {
#endif
  
	vcn_model_t* vcn_model_get_combination(const vcn_model_t *const model1,
					       const vcn_model_t *const model2,
					       double min_length_x_segment);

	vcn_model_t* vcn_model_get_intersection(const vcn_model_t *const model1,
						const vcn_model_t *const model2,
						double min_length_x_segment);

	vcn_model_t* vcn_model_get_union(const vcn_model_t *const model1,
					 const vcn_model_t *const model2,
					 double min_length_x_segment);

#ifdef __cplusplus
}
#endif

#endif
