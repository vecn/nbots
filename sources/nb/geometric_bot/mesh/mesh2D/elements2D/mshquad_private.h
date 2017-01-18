#ifndef _NB_GEOMETRIC_BOT_MESH_MESH2D_ELEMENTS2D_MSHQUAD_PRIVATE_H_
#define _NB_GEOMETRIC_BOT_MESH_MESH2D_ELEMENTS2D_MSHQUAD_PRIVATE_H_

#include <stdint.h>

void nb_mshquad_node_move_x(void *msh, uint32_t i, double x);
void nb_mshquad_node_move_y(void *msh, uint32_t i, double y);
void nb_mshquad_elem_move_x(void *msh, uint32_t i, double x);
void nb_mshquad_elem_move_y(void *msh, uint32_t i, double y);
uint32_t nb_mshquad_get_N_total_adj(const void *msh);

#endif
