#ifndef __NB_GEOMETRIC_BOT_MESH_MODULES2D_IMAGE_DENSITY_H__
#define __NB_GEOMETRIC_BOT_MESH_MODULES2D_IMAGE_DENSITY_H__

#include "nb/image_bot.h"

void nb_mesh_set_img_density(nb_mesh_t *mesh,
			      const nb_image_t *const img,
			      double density_volume);
void nb_mesh_clear_img_density(nb_mesh_t *mesh);

#endif
