#ifndef __NB_PDE_BOT_CONTROL_VOLUME_DRAW_H__
#define __NB_PDE_BOT_CONTROL_VOLUME_DRAW_H__

#include "nb/geometric_bot.h"

void nb_cvfa_draw_integration_mesh(const nb_mesh2D_t *const part,
				   const char *filename, int w, int h);

#endif
