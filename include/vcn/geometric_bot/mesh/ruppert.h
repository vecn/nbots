#ifndef __VCN_GEOMETRIC_BOT_MESH_RUPPERT_H__
#define __VCN_GEOMETRIC_BOT_MESH_RUPPERT_H__

#include <stdbool.h>

#include "vcn/geometric_bot/mesh/mesh2D.h"

void vcn_ruppert_refine(vcn_mesh_t *mesh);
bool vcn_ruppert_insert_vtx(vcn_mesh_t *restrict mesh, const double vertex[2]);

#endif
