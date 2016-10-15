#ifndef __NB_GEOMETRIC_BOT_MESH_CONSTRAINED_DELAUNAY_H__
#define __NB_GEOMETRIC_BOT_MESH_CONSTRAINED_DELAUNAY_H__

#include "nb/geometric_bot/mesh/mesh2D.h"

void nb_tessellator2D_get_constrained_delaunay(nb_tessellator2D_t *mesh,
				      uint32_t N_vertices,
				      const double *const vertices,
				      uint32_t N_segments,
				      const uint32_t *const segments);

#endif
