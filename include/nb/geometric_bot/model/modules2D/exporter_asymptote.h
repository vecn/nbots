#ifndef __NB_GEOMETRIC_BOT_MODEL_MODULES2D_EXPORTER_ASYMPTOTE_H__
#define __NB_GEOMETRIC_BOT_MODEL_MODULES2D_EXPORTER_ASYMPTOTE_H__

#include <stdbool.h>
#include <stdint.h>

#include "nb/geometric_bot/model/model2D.h"

  
/**
 * @brief Export the <a href="http://asymptote.sourceforge.net/">Asymptote</a>
 * code to generate an image of the model.
 * @param[in] model Model to be displayed in the image.
 * @param[in] filename Name of the ASY file.
 */
void vcn_model_export_to_asymptote(const vcn_model_t *const model,
				   const char* filename);

#endif
