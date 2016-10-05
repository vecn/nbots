#ifndef __NB_GEOMETRIC_BOT_MESH_DEWALL_H__
#define __NB_GEOMETRIC_BOT_MESH_DEWALL_H__

#include "nb/geometric_bot/mesh/mesh2D.h"


/**
 * @brief Create a Delaunay triangulation from a given set of vertices.
 * @param[in,out] Mesh structure to store Delaunay triangulation.
 * @param[in] N_vertices Number of vertices.
 * @param[in] vertices Concatenated coordinates of the vertices.
 */
void nb_mesh_get_delaunay(nb_mesh_t *mesh, uint32_t N_vertices,
			  const double *const vertices);

#endif
