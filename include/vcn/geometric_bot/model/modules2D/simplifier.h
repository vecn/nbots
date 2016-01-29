#ifndef __VCN_MODEL2D_SIMPLIFIER_H__
#define __VCN_MODEL2D_SIMPLIFIER_H__

#include <stdint.h>
#include "vcn/container_bot/container.h"
#include "vcn/geometric_bot/model/model2D.h"

#ifdef __cplusplus
extern "C" {
#endif
	vcn_container_t* vcn_model_generate_wires
					(const vcn_model_t *const model);
	void vcn_model_collapse_small_segments(vcn_model_t* model,
					       double tolerance,
					       uint32_t N_fixed_vertices,
					       uint32_t* fixed_vertices);

	void vcn_model_collapse_colinear_vertices(vcn_model_t* model,
						  uint32_t N_fixed_vertices,
						  uint32_t* fixed_vertices,
						  double tolerance);

	void vcn_model_unify_edge(vcn_model_t* model, double* vtx1, double* vtx2);

#ifdef __cplusplus
}
#endif

#endif
