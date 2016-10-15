#ifndef __NB_GEOMETRIC_BOT_MESH_MODULES2D_IMAGE_DENSITY_H__
#define __NB_GEOMETRIC_BOT_MESH_MODULES2D_IMAGE_DENSITY_H__

#include "nb/image_bot.h"

void nb_tessellator2D__set_img_density(nb_tessellator2D__t *mesh,
			      const nb_image_t *const img,
			      double density_volume);
void nb_tessellator2D__clear_img_density(nb_tessellator2D__t *mesh);

#endif
