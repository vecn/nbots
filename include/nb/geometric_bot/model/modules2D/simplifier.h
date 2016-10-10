#ifndef __NB_GEOMETRIC_BOT_MODEL_MODULES2D_SIMPLIFIER_H__
#define __NB_GEOMETRIC_BOT_MODEL_MODULES2D_SIMPLIFIER_H__

#include <stdint.h>
#include "nb/container_bot/container.h"
#include "nb/geometric_bot/model/model2D.h"

nb_container_t* nb_model_generate_wires(const nb_model_t *const model);

void nb_model_collapse_small_segments(nb_model_t* model,
				       double tolerance,
				       uint32_t N_fixed_vertices,
				       uint32_t* fixed_vertices);

void nb_model_collapse_colinear_vertices(nb_model_t* model,
					  uint32_t N_fixed_vertices,
					  uint32_t* fixed_vertices,
					  double tolerance);

void nb_model_unify_edge(nb_model_t* model, double* vtx1, double* vtx2);

#endif
