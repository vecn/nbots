#ifndef __VCN_GEOMETRIC_BOT_MESH_DEWALL_H__
#define __VCN_GEOMETRIC_BOT_MESH_DEWALL_H__

#include "vcn/geometric_bot/mesh/mesh2D.h"

#ifdef __cplusplus
extern "C" {
#endif

	/**
	 * @brief Create a Delaunay triangulation from a given set of vertices.
	 * @param[in] N_vertices Number of vertices.
	 * @param[in] vertices Concatenated coordinates of the vertices.
	 * @return The mesh corresponding to the triangulation if success, or
	 * NULL if something goes wrong.
	 */
	vcn_mesh_t* vcn_dewall_get_delaunay(uint32_t N_vertices,
					    const double *const vertices);
#ifdef __cplusplus
}
#endif

#endif
