#ifndef __NB_GEOMETRIC_BOT_MESH2D_FILE_FORMAT_NBT_H__
#define __NB_GEOMETRIC_BOT_MESH2D_FILE_FORMAT_NBT_H__

#include "nb/geometric_bot.h"

int nb_mesh2D_save_nbt(const nb_mesh2D_t *mesh,  const char *name);
int nb_mesh2D_read_type_nbt(const char *name, nb_mesh2D_type *type);
int nb_mesh2D_read_nbt(nb_mesh2D_t *mesh, const char *name);

#endif
