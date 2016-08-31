#ifndef __NB_GEOMETRIC_BOT_MESH_RUPPERT_H__
#define __NB_GEOMETRIC_BOT_MESH_RUPPERT_H__

#include <stdbool.h>

#include "nb/geometric_bot/mesh/mesh2D.h"

void vcn_ruppert_refine(vcn_mesh_t *mesh);
bool vcn_ruppert_insert_vtx(vcn_mesh_t *restrict mesh,
			    const double vertex[2]);

void nb_ruppert_split_all_subsgm(nb_mesh_t *mesh);
void nb_ruppert_split_trg_with_all_nodes_in_sgm(nb_mesh_t *mesh);

#endif
