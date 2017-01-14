#ifndef _NB_GEOMETRIC_BOT_MESH_MESH2D_ELEMENTS2D_MSHPACK_PRIVATE_H_
#define _NB_GEOMETRIC_BOT_MESH_MESH2D_ELEMENTS2D_MSHPACK_PRIVATE_H_

#include <stdint.h>

void nb_mshpack_node_move_x(void *msh, uint32_t i, double x);
void nb_mshpack_node_move_y(void *msh, uint32_t i, double y);
void nb_mshpack_elem_move_x(void *msh, uint32_t i, double x);
void nb_mshpack_elem_move_y(void *msh, uint32_t i, double y);
uint32_t nb_mshpack_get_N_total_adj(const void *msh);

#endif
