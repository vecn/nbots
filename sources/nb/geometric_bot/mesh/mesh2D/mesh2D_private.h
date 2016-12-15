#ifndef __NB_GEOMETRIC_BOT_MESH_MESH2D_PRIVATE_H__
#define __NB_GEOMETRIC_BOT_MESH_MESH2D_PRIVATE_H__

#include <stdint.h>

#include "mesh2D_struct.h"

typedef struct {
	void (*node_move_x)(void *msh, uint32_t i, double u);
	void (*node_move_y)(void *msh, uint32_t i, double v);
	void (*elem_move_x)(void *msh, uint32_t i, double u);
	void (*elem_move_y)(void *msh, uint32_t i, double v);
} nb_mesh2D_private_i;

void nb_mesh2D_init_private_interface(nb_mesh2D_private_i *priv,
					 const nb_mesh2D_t *part);

#endif
