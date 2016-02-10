#ifndef __NB_GEOMETRIC_BOT_MESH_CONSTRAINED_DELAUNAY_H__
#define __NB_GEOMETRIC_BOT_MESH_CONSTRAINED_DELAUNAY_H__

#include "nb/geometric_bot/mesh/mesh2D.h"

void vcn_mesh_get_constrained_delaunay(vcn_mesh_t *mesh,
				       uint32_t N_vertices,
				       const double *const vertices,
				       uint32_t N_segments,
				       const uint32_t *const segments);

#endif
