#ifndef __NB_GEOMETRIC_BOT_MESH2D_FILE_FORMAT_VTK_H__
#define __NB_GEOMETRIC_BOT_MESH2D_FILE_FORMAT_VTK_H__

#include "nb/geometric_bot.h"

int nb_mesh2D_save_vtk(const nb_mesh2D_t *mesh,  const char *name);
int nb_mesh2D_read_vtk(nb_mesh2D_t *mesh, const char *name);

#endif
