#ifndef __NB_GEOMETRIC_BOT_MESH_PARTITION_H__
#define __NB_GEOMETRIC_BOT_MESH_PARTITION_H__

#include <stdint.h>

#include "nb/graph_bot.h"
#include "nb/geometric_bot/mesh/mesh2D.h"

typedef enum {
	NB_TRIAN, NB_QUAD, NB_POLY, NB_DISK
} nb_partition_type;

typedef struct nb_partition_s nb_partition_t;

uint32_t nb_partition_get_memsize(nb_partition_type  type);
void nb_partition_init(nb_partition_t *part, nb_partition_type  type);
void nb_partition_copy(nb_partition_t *part, const nb_partition_t* srcpart);
void nb_partition_finish(nb_partition_t *part);
nb_partition_t* nb_partition_create(nb_partition_type type);
nb_partition_t* nb_partition_clone(nb_partition_t* part);
void nb_partition_clear(nb_partition_t* part);
void nb_partition_destroy(nb_partition_t* part);

nb_partition_type nb_partition_get_type(const nb_partition_t *part);

uint32_t nb_partition_get_N_invtx(const nb_partition_t *part);
uint32_t nb_partition_get_N_insgm(const nb_partition_t *part);
uint32_t nb_partition_get_N_nodes(const nb_partition_t *part);
uint32_t nb_partition_get_N_edges(const nb_partition_t *part);
uint32_t nb_partition_get_N_elems(const nb_partition_t *part);

double nb_partition_node_get_x(const nb_partition_t *part, uint32_t id);
double nb_partition_node_get_y(const nb_partition_t *part, uint32_t id);

uint32_t nb_partition_edge_get_1n(const nb_partition_t *part, uint32_t id);
uint32_t nb_partition_edge_get_2n(const nb_partition_t *part, uint32_t id);

double nb_partition_elem_get_x(const nb_partition_t *part, uint32_t id);
double nb_partition_elem_get_y(const nb_partition_t *part, uint32_t id);
double nb_partition_elem_get_area(const nb_partition_t *part, uint32_t id);
double nb_partition_elem_face_get_length(const nb_partition_t *part,
					 uint32_t elem_id, uint16_t face_id);

uint32_t nb_partition_elem_get_N_adj(const nb_partition_t *part, uint32_t id);
uint32_t nb_partition_elem_get_adj(const nb_partition_t *part,
				   uint32_t elem_id, uint8_t adj_id);
uint32_t nb_partition_elem_get_N_ngb(const nb_partition_t *part, uint32_t id);
uint32_t nb_partition_elem_get_ngb(const nb_partition_t *part,
				   uint32_t elem_id, uint8_t ngb_id);
bool nb_partition_elem_has_ngb(const nb_partition_t *part, uint32_t elem_id,
			       uint16_t face_id);
uint32_t nb_partition_get_invtx(const nb_partition_t *part, uint32_t id);
uint32_t nb_partition_insgm_get_N_subsgm(const nb_partition_t *part,
					 uint32_t id);
uint32_t nb_partition_get_N_nodes_x_insgm(const nb_partition_t *part,
					  uint32_t id);
uint32_t nb_partition_get_node_x_insgm(const nb_partition_t *part,
				       uint32_t sgm_id, uint32_t node_id);
double nb_partition_insgm_get_length(const nb_partition_t *part,
				     uint32_t sgm_id);
double nb_partition_insgm_subsgm_get_length(const nb_partition_t *part,
					    uint32_t sgm_id,
					    uint32_t subsgm_id);
void nb_partition_insgm_get_elem_adj(const nb_partition_t *part,
				     uint32_t **elem_adj);

void nb_partition_load_elem_graph(const nb_partition_t *part,
				  vcn_graph_t *graph);
void nb_partition_load_nodal_graph(const nb_partition_t *part,
				   vcn_graph_t *graph);
void nb_partition_load_interelem_graph(const nb_partition_t *part,
				       vcn_graph_t *graph);
void nb_partition_load_from_mesh(nb_partition_t *part,
				 nb_mesh_t *mesh);
void nb_partition_get_enveloping_box(const nb_partition_t *part,
				     double box[4]);
bool nb_partition_is_vtx_inside(const nb_partition_t *part,
				double x, double y);
void nb_partition_draw(const nb_partition_t *part, const char *filename,
		       int width, int height);
void nb_partition_build_model(const nb_partition_t *part, nb_model_t *model);

/**
 * @brief Build a geometry model from a triangulation with disabled
 * triangles.
 * @param[in] part Input triangulation.
 * @param[in] elems_enabled Array of booleans indicating wich triangles are
 * enabled. NULL to enable all.
 * @param[out] model Model generated.
 * @param[out] N_input_vtx Stores the number of vertices forming
 * the boundary that appears in the original input.
 * @param[out] input_vtx Stores an array of the vertices ids in
 * the new model corresponding to original input vertices.
 */
void nb_partition_build_model_disabled_elems
			(const nb_partition_t *part,
			   const bool *elems_enabled,
			   nb_model_t *model,
			   uint32_t *N_input_vtx,
			   uint32_t **input_vtx);


#endif
