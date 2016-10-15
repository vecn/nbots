#ifndef __NB_GEOMETRIC_BOT_MESH_RUPPERT_H__
#define __NB_GEOMETRIC_BOT_MESH_RUPPERT_H__

#include <stdbool.h>

#include "nb/geometric_bot/mesh/mesh2D.h"

#include "../mesh2D_structs.h"

void nb_ruppert_refine(nb_tessellator2D_t *mesh);
bool nb_ruppert_insert_vtx(nb_tessellator2D_t *restrict mesh,
			    const double vertex[2]);
void nb_ruppert_insert_verified_vtx(nb_tessellator2D_t *restrict mesh,
				    msh_trg_t *trg,
				    const msh_vtx_t *vtx);
void nb_ruppert_insert_verified_subsgm_midpoint(nb_tessellator2D_t *restrict mesh,
						msh_edge_t *sgm);

#endif
