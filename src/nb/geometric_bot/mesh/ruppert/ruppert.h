#ifndef __NB_GEOMETRIC_BOT_MESH_RUPPERT_H__
#define __NB_GEOMETRIC_BOT_MESH_RUPPERT_H__

#include <stdbool.h>

#include "nb/geometric_bot/mesh/mesh2D.h"

#include "../mesh2D_structs.h"

void vcn_ruppert_refine(vcn_mesh_t *mesh);
bool vcn_ruppert_insert_vtx(vcn_mesh_t *restrict mesh,
			    const double vertex[2]);
void nb_ruppert_insert_verified_vtx(nb_mesh_t *restrict mesh,
				    msh_trg_t *trg,
				    const msh_vtx_t *vtx);
void nb_ruppert_insert_verified_subsgm_midpoint(nb_mesh_t *restrict mesh,
						msh_edge_t *sgm);

#endif
