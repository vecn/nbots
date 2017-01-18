#ifndef _NB_GEOMETRIC_BOT_MESH_MESH2D_ELEMENTS2D_MSHQUAD_FILE_FORMAT_NBT_H_
#define _NB_GEOMETRIC_BOT_MESH_MESH2D_ELEMENTS2D_MSHQUAD_FILE_FORMAT_NBT_H_

#include <stdio.h>

#include "nb/io_bot.h"

void nb_mshquad_write_data_nbt(FILE *fp, const void *msh);
int nb_mshquad_read_data_nbt(nb_cfreader_t *cfr, void *msh);

#endif
