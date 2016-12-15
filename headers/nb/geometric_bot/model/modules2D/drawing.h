#ifndef __NB_GEOMETRIC_BOT_MODEL_MODULES2D_DRAWING_H__
#define __NB_GEOMETRIC_BOT_MODEL_MODULES2D_DRAWING_H__

#include <stdint.h>
#include <stdbool.h>
#include "nb/geometric_bot.h"

void nb_model_draw(const nb_model_t *const model,
		    const char* filename,
		    int width, int height);

#endif
