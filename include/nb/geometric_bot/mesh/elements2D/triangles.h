#ifndef __NB_GEOMETRIC_BOT_MESH_ELEMENTS2D_TRIANGLES_H__
#define __NB_GEOMETRIC_BOT_MESH_ELEMENTS2D_TRIANGLES_H__

#include <stdbool.h>
#include <stdint.h>
#include "triangles_struct.h"
#include "nb/geometric_bot/mesh/mesh2D.h"
#include "nb/graph_bot.h"

uint32_t vcn_msh3trg_get_memsize(void);
void vcn_msh3trg_init(vcn_msh3trg_t *msh3trg);
void vcn_msh3trg_finish(vcn_msh3trg_t *msh3trg);

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

bool nb_msh3trg_is_vtx_inside(const vcn_msh3trg_t *msh3trg,
			      const double x[2]);

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


void vcn_msh3trg_load_from_mesh(vcn_msh3trg_t *msh3trg,
				const vcn_mesh_t *const mesh);

void vcn_msh3trg_relabel(vcn_msh3trg_t* msh3trg,
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
