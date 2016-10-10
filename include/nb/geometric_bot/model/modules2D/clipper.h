#ifndef __NB_GEOMETRIC_BOT_MODEL_MODULES2D_CLIPPER_H__
#define __NB_GEOMETRIC_BOT_MODEL_MODULES2D_CLIPPER_H__

#include <stdbool.h>
#include <stdint.h>

#include "nb/geometric_bot/model/model2D.h"
#include "nb/graph_bot.h"

void nb_model_get_combination(nb_model_t *model,
			       const nb_model_t *const model1,
			       const nb_model_t *const model2);

void nb_model_get_intersection(nb_model_t *model,
				const nb_model_t *const model1,
				const nb_model_t *const model2);

void nb_model_get_union(nb_model_t *model,
			 const nb_model_t *const model1,
			 const nb_model_t *const model2);

void nb_model_get_difference(nb_model_t *model,
			      const nb_model_t *const model1,
			      const nb_model_t *const model2);

void nb_model_get_substraction(nb_model_t *model,
				const nb_model_t *const model1,
				const nb_model_t *const model2);

#endif
