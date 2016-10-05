#ifndef __NB_GEOMETRIC_BOT_MESH_PARTITION_ELEMENTS2D_MSH3TRG_H__
#define __NB_GEOMETRIC_BOT_MESH_PARTITION_ELEMENTS2D_MSH3TRG_H__

#include <stdbool.h>
#include <stdint.h>

#include "nb/geometric_bot/model/model2D.h"
#include "nb/geometric_bot/mesh/mesh2D.h"
#include "nb/geometric_bot/mesh/partition/info.h"

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
void nb_msh3trg_edge_get_midpoint(const void *msh,
				  uint32_t face_id, double w,
				  double midpoint[2]);
double nb_msh3trg_edge_get_normal(const void *msh, uint32_t face_id,
				  double normal[2]);
double nb_msh3trg_elem_get_x(const void *msh, uint32_t id);
double nb_msh3trg_elem_get_y(const void *msh, uint32_t id);
double nb_msh3trg_elem_get_area(const void *msh, uint32_t id);
double nb_msh3trg_elem_get_radius(const void *msh, uint32_t id);
double nb_msh3trg_elem_get_apotem(const void *msh, uint32_t id);
uint32_t nb_msh3trg_elem_find_edge(const void *msh, uint32_t id,
				   uint16_t local_face_id);
double nb_msh3trg_elem_face_get_length(const void *msh,
				       uint32_t elem_id,
				       uint16_t face_id);
double nb_msh3trg_elem_face_get_normal(const void *msh, uint32_t elem_id,
				       uint16_t face_id, double normal[2]);
double nb_msh3trg_elem_ngb_get_normal(const void *msh,
				      uint32_t elem_id, uint16_t face_id,
				      double normal[2]);
uint32_t nb_msh3trg_elem_get_N_adj(const void *msh, uint32_t id);
uint32_t nb_msh3trg_elem_get_adj(const void *msh,
				 uint32_t elem_id, uint8_t adj_id);
uint32_t nb_msh3trg_elem_get_ngb(const void *msh,
				 uint32_t elem_id, uint8_t ngb_id);
bool nb_msh3trg_elem_has_ngb(const void *msh, uint32_t elem_id,
			     uint16_t ngb_id);
uint32_t nb_msh3trg_get_invtx(const void *msh, uint32_t id);
uint32_t nb_msh3trg_insgm_get_N_nodes(const void *msh, uint32_t id);
uint32_t nb_msh3trg_insgm_get_node(const void *msh, uint32_t sgm_id,
				     uint32_t node_id);

void nb_msh3trg_load_from_mesh(void *msh3trg, nb_mesh_t *mesh);
void nb_msh3trg_set_nodal_permutation(void *msh, const uint32_t *perm);
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

void nb_msh3trg_build_model(const void *msh3trg, nb_model_t *model);

void nb_msh3trg_build_model_disabled_elems(const void *msh3trg_ptr,
					   const bool *elems_enabled,
					   nb_model_t *model,
					   uint32_t *N_input_vtx,
					   uint32_t **input_vtx);

#endif
