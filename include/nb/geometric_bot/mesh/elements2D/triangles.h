#ifndef __NB_GEOMETRIC_BOT_MESH_ELEMENTS2D_TRIANGLES_H__
#define __NB_GEOMETRIC_BOT_MESH_ELEMENTS2D_TRIANGLES_H__

#include <stdbool.h>
#include <stdint.h>
#include "nb/geometric_bot/mesh/mesh2D.h"
#include "nb/graph_bot.h"

uint32_t nb_msh3trg_get_memsize(void);
void nb_msh3trg_init(void *msh3trg);
void nb_msh3trg_finish(void *msh3trg);
void nb_msh3trg_copy(void *msh3trg, const void *src);

void* nb_msh3trg_create(void);
void* nb_msh3trg_clone(void* msh3trg);
void nb_msh3trg_clear(void* msh3trg);
void nb_msh3trg_destroy(void* msh3trg);

uint32_t nb_msh3trg_get_N_invtx(const void *msh);
uint32_t nb_msh3trg_get_N_insgm(const void *msh);
uint32_t nb_msh3trg_get_N_nodes(const void *msh);
uint32_t nb_msh3trg_get_N_edges(const void *msh);
uint32_t nb_msh3trg_get_N_elems(const void *msh);
double nb_msh3trg_get_x_node(const void *msh, uint32_t id);
double nb_msh3trg_get_y_node(const void *msh, uint32_t id);
uint32_t nb_msh3trg_get_1n_edge(const void *msh, uint32_t id);
uint32_t nb_msh3trg_get_2n_edge(const void *msh, uint32_t id);
double nb_msh3trg_get_x_elem(const void *msh, uint32_t id);
double nb_msh3trg_get_y_elem(const void *msh, uint32_t id);
uint32_t nb_msh3trg_elem_get_N_adj(const void *msh, uint32_t id);
uint32_t nb_msh3trg_elem_get_adj(const void *msh,
				 uint32_t elem_id, uint8_t adj_id);
uint32_t nb_msh3trg_elem_get_N_ngb(const void *msh, uint32_t id);
uint32_t nb_msh3trg_elem_get_ngb(const void *msh,
				 uint32_t elem_id, uint8_t ngb_id);
uint32_t nb_msh3trg_get_invtx(const void *msh, uint32_t id);
uint32_t nb_msh3trg_get_N_nodes_x_insgm(const void *msh, uint32_t id);
uint32_t nb_msh3trg_get_node_x_insgm(const void *msh, uint32_t sgm_id,
				     uint32_t node_id);

bool nb_msh3trg_is_vtx_inside(const void *msh3trg,
			      const double x[2]);

/**
 * @brief Build a graph from the vertices of the mesh.
 * @param[in] msh3trg Mesh representing the graph.
 */
void nb_msh3trg_load_vtx_graph(const void *const msh3trg,
			       nb_graph_t *graph);

/**
 * @brief Build a graph from the elements mesh.
 * @param[in] msh3trg Mesh representing the graph.
 */
void nb_msh3trg_load_elem_graph(const void *const msh3trg,
			       nb_graph_t *graph);


void nb_msh3trg_load_from_mesh(void *msh3trg,
				const vcn_mesh_t *const mesh);

void nb_msh3trg_relabel(void* msh3trg,
			 uint32_t* (*labeling)(const nb_graph_t *const));

/**
 * @brief Disable the elements sharing single point connections.
 * @param[in] msh3trg Triangular mesh, array with indicators of enabled and
 * disabled triangular elements.
 * @param[out] enabled_elements Turn to false the indices of this array
 * corresponding to elements which share single point connections.
 */
void nb_msh3trg_disable_single_point_connections
(const void *const msh3trg, bool* enabled_elements);

void nb_msh3trg_get_enveloping_box(const void *msh3trg_ptr);

#endif
