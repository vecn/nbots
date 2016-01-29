#ifndef __VCN_MODEL2D_EXPORTER_ASYMPTOTE_H__
#define __VCN_MODEL2D_EXPORTER_ASYMPTOTE_H__

#include <stdbool.h>
#include <stdint.h>

#include "vcn/geometric_bot/model/model2D.h"

#ifdef __cplusplus
extern "C" {
#endif
  
	/**
	 * @brief Export the <a href="http://asymptote.sourceforge.net/">Asymptote</a>
	 * code to generate an image of the model.
	 * @param[in] model Model to be displayed in the image.
	 * @param[in] filename Name of the ASY file.
	 */
	void vcn_model_export_to_asymptote(const vcn_model_t *const model,
					   const char* filename);

#ifdef __cplusplus
}
#endif

#endif
