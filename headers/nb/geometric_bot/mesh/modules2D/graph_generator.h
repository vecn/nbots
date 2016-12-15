#ifndef __NB_GEOMETRIC_BOT_MESH_MODULES2D_GRAPH_GENERATOR_H__
#define __NB_GEOMETRIC_BOT_MESH_MODULES2D_GRAPH_GENERATOR_H__

#include "nb/graph_bot.h"
#include "nb/geometric_bot/mesh/tessellator2D.h"

/**
 * @brief Build a graph from the triangular mesh. Each element of the mesh
 * represent a node in the graph.
 * @n<b>WARNING:</b> Assuming that triangles has an ID as attribute.
 * @param[in] mesh Mesh representing the graph.
 * @return The graph generated from the elements (triangles) of the mesh.
 */
void nb_tessellator2D_load_elem_graph(const nb_tessellator2D_t *const mesh,
				      nb_graph_t *graph);

#endif
