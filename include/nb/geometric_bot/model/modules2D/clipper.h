#ifndef __NB_GEOMETRIC_BOT_MODEL_MODULES2D_CLIPPER_H__
#define __NB_GEOMETRIC_BOT_MODEL_MODULES2D_CLIPPER_H__

#include <stdbool.h>
#include <stdint.h>

#include "nb/geometric_bot/mesh/elements2D/triangles_struct.h"
#include "nb/geometric_bot/model/model2D.h"
#include "nb/graph_bot.h"

void vcn_model_get_combination(vcn_model_t *model,
			       const vcn_model_t *const model1,
			       const vcn_model_t *const model2);

void vcn_model_get_intersection(vcn_model_t *model,
				const vcn_model_t *const model1,
				const vcn_model_t *const model2);

void vcn_model_get_union(vcn_model_t *model,
			 const vcn_model_t *const model1,
			 const vcn_model_t *const model2);

void vcn_model_get_difference(vcn_model_t *model,
			      const vcn_model_t *const model1,
			      const vcn_model_t *const model2);

void vcn_model_get_substraction(vcn_model_t *model,
				const vcn_model_t *const model1,
				const vcn_model_t *const model2);

#endif
