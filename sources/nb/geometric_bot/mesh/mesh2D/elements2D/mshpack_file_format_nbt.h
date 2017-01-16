#ifndef _NB_GEOMETRIC_BOT_MESH_MESH2D_ELEMENTS2D_MSHPACK_FILE_FORMAT_NBT_H_
#define _NB_GEOMETRIC_BOT_MESH_MESH2D_ELEMENTS2D_MSHPACK_FILE_FORMAT_NBT_H_

#include <stdio.h>

#include "nb/cfreader_bot.h"

void nb_mshpack_write_data_nbt(FILE *fp, const void *msh);
int nb_mshpack_read_data_nbt(nb_cfreader_t *cfr, void *msh);

#endif
