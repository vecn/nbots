#ifndef __NB_GEOMETRIC_BOT_MESH_MODULES2D_GRAPH_GENERATOR_H__
#define __NB_GEOMETRIC_BOT_MESH_MODULES2D_GRAPH_GENERATOR_H__

/**
 * @brief Build a graph from the triangular mesh. Each vertex of the mesh
 * represent a node in the graph.
 * @n<b>WARNING:</b> Assuming that vertices has an ID as attribute.
 * @param[in] mesh Mesh representing the graph.
 * @return The graph generated from the vertices of the mesh.
 */
nb_graph_t* nb_mesh_create_vtx_graph(const nb_mesh_t *const mesh);

/**
 * @brief Build a graph from the triangular mesh. Each element of the mesh
 * represent a node in the graph.
 * @n<b>WARNING:</b> Assuming that triangles has an ID as attribute.
 * @param[in] mesh Mesh representing the graph.
 * @return The graph generated from the elements (triangles) of the mesh.
 */
nb_graph_t* nb_mesh_create_elem_graph(const nb_mesh_t *const mesh);

#endif
