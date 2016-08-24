#ifndef __NB_GEOMETRIC_BOT_MESH_ELEMENTS2D_TRIANGLES_H__
#define __NB_GEOMETRIC_BOT_MESH_ELEMENTS2D_TRIANGLES_H__

#include <stdbool.h>
#include <stdint.h>
#include "nb/geometric_bot/mesh/mesh2D.h"
#include "nb/graph_bot.h"

typedef struct nb_msh3trg_s nb_msh3trg_t;

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
double nb_msh3trg_node_get_x(const void *msh, uint32_t id);
double nb_msh3trg_node_get_y(const void *msh, uint32_t id);
uint32_t nb_msh3trg_edge_get_1n(const void *msh, uint32_t id);
uint32_t nb_msh3trg_edge_get_2n(const void *msh, uint32_t id);
double nb_msh3trg_elem_get_x(const void *msh, uint32_t id);
double nb_msh3trg_elem_get_y(const void *msh, uint32_t id);
double nb_msh3trg_elem_get_area(const void *msh, uint32_t id);
double nb_msh3trg_elem_face_get_length(const void *msh,
				       uint32_t elem_id,
				       uint16_t face_id);
uint32_t nb_msh3trg_elem_get_N_adj(const void *msh, uint32_t id);
uint32_t nb_msh3trg_elem_get_adj(const void *msh,
				 uint32_t elem_id, uint8_t adj_id);
uint32_t nb_msh3trg_elem_get_N_ngb(const void *msh, uint32_t id);
uint32_t nb_msh3trg_elem_get_ngb(const void *msh,
				 uint32_t elem_id, uint8_t ngb_id);
bool nb_msh3trg_elem_has_ngb(const void *msh, uint32_t elem_id,
			     uint16_t ngb_id);
uint32_t nb_msh3trg_get_invtx(const void *msh, uint32_t id);
uint32_t nb_msh3trg_get_N_nodes_x_insgm(const void *msh, uint32_t id);
uint32_t nb_msh3trg_get_node_x_insgm(const void *msh, uint32_t sgm_id,
				     uint32_t node_id);

void nb_msh3trg_load_elem_graph(const void *const msh3trg,
				nb_graph_t *graph);

void nb_msh3trg_load_nodal_graph(const void *const msh3trg,
				 nb_graph_t *graph);

void nb_msh3trg_load_interelem_graph(const void *const msh3trg,
				     nb_graph_t *graph);


void nb_msh3trg_load_from_mesh(void *msh3trg, vcn_mesh_t *mesh);

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

void nb_msh3trg_get_enveloping_box(const void *msh3trg_ptr, double box[4]);

bool nb_msh3trg_is_vtx_inside(const void *msh3trg, double x, double y);

void nb_msh3trg_draw(const void *msh3trg_ptr,
		     const char* filename, int width, int height);

/**
 * @brief Export a PNG image of the partition ans its subdomains.
 * @param[in] part Partition to be displayed.
 * @param[in] k_part Number of partitions.
 * @param[in] part Partition corresponding to the vertex.
 * @param[in] filename Name of the image file.
 * @param[in] width Image width.
 * @param[in] height Image height.
 * @param[in] k_to_draw Partition to draw. Zero to draw all of them.
 * @param[in] scale_partitions Scale of the size of the partitions, in (0,1].
 */
void nb_msh3trg_draw_subdomain(const void *msh3trg_ptr,
			       const char* filename, int width, int height,
			       uint32_t k_part, const uint32_t *const part,
			       uint32_t k_to_draw, double scale_subdomains);

void nb_msh3trg_build_model(const void *msh3trg, nb_model_t *model);

void nb_msh3trg_build_model_disabled_elems(const void *msh3trg_ptr,
					   const bool *elems_enabled,
					   nb_model_t *model,
					   uint32_t *N_input_vtx,
					   uint32_t **input_vtx);
#endif
