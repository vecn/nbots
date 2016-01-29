#ifndef __VCN_GEOMETRIC_BOT_MESH_CONSTRAINED_DELAUNAY_H__
#define __VCN_GEOMETRIC_BOT_MESH_CONSTRAINED_DELAUNAY_H__

#include "vcn/geometric_bot/mesh/mesh2D.h"

#ifdef __cplusplus
extern "C" {
#endif

	vcn_mesh_t* vcn_mesh_get_constrained_delaunay
					(uint32_t N_vertices,
					 const double *const vertices,
					 uint32_t N_segments,
					 const uint32_t *const segments);

#ifdef __cplusplus
}
#endif

#endif
