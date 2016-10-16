#ifndef __NB_GEOMETRIC_BOT_MESH_MODULES2D_IMAGE_DENSITY_H__
#define __NB_GEOMETRIC_BOT_MESH_MODULES2D_IMAGE_DENSITY_H__

#include "nb/image_bot.h"

#include "nb/geometric_bot/mesh/tessellator2D.h"

void nb_tessellator2D_set_img_density(nb_tessellator2D_t *mesh,
				      const nb_image_t *const img,
				      double density_volume);
void nb_tessellator2D_clear_img_density(nb_tessellator2D_t *mesh);

#endif
