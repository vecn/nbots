#ifndef __VCN_GEOMETRIC_BOT_MESH_ELEMENTS2D_TRIANGLES_H__
#define __VCN_GEOMETRIC_BOT_MESH_ELEMENTS2D_TRIANGLES_H__

#include <stdbool.h>
#include <stdint.h>
#include "triangles_struct.h"
#include "nb/geometric_bot/mesh/mesh2D.h"
#include "nb/graph_bot.h"

/**
 * @brief Create an empty structure to store triangulations.
 * @return The empty structure.
 */
vcn_msh3trg_t* vcn_msh3trg_create(void);

/**
 * @brief Create an identical copy of the triangulation.
 * @param[in] msh3trg Triangulation to be cloned.
 * @return Cloned triangulation.
 */
vcn_msh3trg_t* vcn_msh3trg_clone(vcn_msh3trg_t* msh3trg);
  
/**
 * @brief Clear the triangulation structure.
 * @param[in] msh3trg Structure to be cleared.
 */
void vcn_msh3trg_clear(vcn_msh3trg_t* msh3trg);

/**
 * @brief Destroy the triangulation structure.
 * @param[in] msh3trg Structure to be destroyed
 */
void vcn_msh3trg_destroy(vcn_msh3trg_t* msh3trg);

/**
 * @brief Build a graph from the vertices of the mesh.
 * @param[in] msh3trg Mesh representing the graph.
 * @return The graph generated from the vertices of the mesh.
 */
vcn_graph_t* vcn_msh3trg_create_vtx_graph(const vcn_msh3trg_t *const msh3trg);

/**
 * @brief Build a graph from the elements mesh.
 * @param[in] msh3trg Mesh representing the graph.
 * @return The graph generated from the elements mesh.
 */
vcn_graph_t* vcn_msh3trg_create_elem_graph(const vcn_msh3trg_t *const msh3trg);

/**
 * @brief Create a read-only triangular mesh from the write-only mesh.
 * @param[in] include_triangles Include or exclude the connectivity 
 * matrix from the output.
 * @param[in] include_input_segments Include or exclude the correspondence-
 * table relating input segments ids with mesh edges ids.
 * @param[in] include_input_vertices Include or exclude the correspondence-
 * table relating input vertices ids with mesh vertices ids.
 * @n<b>WARNING</b>: If a vertex was deleted by the function
 * vcn_mesh_delete_isolated_vertices(), the corresponding id will be
 * out of bounds.
 * @param[in] include_neighbours Include or exclude the triangle's neighbours
 * relation from the output.
 * @param[in] labeling Labeling function. This function computes the
 * best permutation of labels for a particular operation, such as LU decomposition
 * or to minimize the band-width of the corresponding sparse-matrix.
 * NULL to generate an arbitrary labeling.
 * @return The triangular mesh if success, NULL if something goes wrong.
 */
vcn_msh3trg_t* vcn_mesh_get_msh3trg(const vcn_mesh_t *const mesh,
				    bool include_edges,
				    bool include_triangles,
				    bool include_neighbours,
				    bool include_input_vertices,
				    bool include_input_segments);

void vcn_mesh_trg2D_relabel(vcn_msh3trg_t* msh3trg,
			    uint32_t* (*labeling)(const vcn_graph_t *const));

/**
 * @brief Disable the elements sharing single point connections.
 * @param[in] msh3trg Triangular mesh, array with indicators of enabled and
 * disabled triangular elements.
 * @param[out] enabled_elements Turn to false the indices of this array
 * corresponding to elements which share single point connections.
 */
void vcn_msh3trg_disable_single_point_connections
(const vcn_msh3trg_t *const msh3trg, bool* enabled_elements);

#endif
