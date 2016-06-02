#ifndef __NB_GEOMETRIC_BOT_MODEL_MODULES2D_DRAWING_H__
#define __NB_GEOMETRIC_BOT_MODEL_MODULES2D_DRAWING_H__

#include <stdint.h>
#include <stdbool.h>
#include "nb/geometric_bot.h"

/**
 * @brief Export a PNG image of the model.
 * @param[in] model Model to be displayed in the image.
 * @param[in] filename Name of the PNG file.
 * @param[in] width Width of the image.
 * @param[in] height of the image.
 * segments.
 */
void vcn_model_export_png(const vcn_model_t *const model,
			  const char* filename,
			  int width, int height);

/**
 * @brief Export an EPS image of the model.
 * @param[in] model Model to be displayed in the image.
 * @param[in] filename Name of the EPS file.
 * @param[in] width Width of the image.
 * @param[in] height of the image.
 * segments.
 */
void vcn_model_export_eps(const vcn_model_t *const model,
			  const char* filename,
			  int width, int height);

#endif
