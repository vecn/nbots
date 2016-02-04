#ifndef __VCN_GEOMETRIC_BOT_MESH_MODULES2D_IMAGE_DENSITY_H__
#define __VCN_GEOMETRIC_BOT_MESH_MODULES2D_IMAGE_DENSITY_H__

#include "vcn/image_bot.h"

void vcn_mesh_set_img_density(vcn_mesh_t *mesh,
			      const vcn_image_t *const img,
			      double scale, double xdisp, double ydisp,
			      double max_density);
void vcn_mesh_clear_img_density(vcn_mesh_t *mesh);

#endif
