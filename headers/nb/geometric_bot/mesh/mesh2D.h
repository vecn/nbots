#ifndef __NB_GEOMETRIC_BOT_MESH2D_H__
#define __NB_GEOMETRIC_BOT_MESH2D_H__

#include <stdint.h>

#include "nb/graph_bot.h"
#include "nb/geometric_bot/mesh/tessellator2D.h"
#include "nb/geometric_bot/mesh/mesh2D/info.h"

typedef enum {
  NB_TRIAN, NB_QUAD, NB_POLY, NB_DISK
} nb_mesh2D_type;

typedef struct nb_mesh2D_s nb_mesh2D_t;

uint32_t nb_mesh2D_get_memsize(nb_mesh2D_type  type);
void nb_mesh2D_init(nb_mesh2D_t *mesh, nb_mesh2D_type  type);
void nb_mesh2D_init_from_msh(nb_mesh2D_t *mesh, void *msh,
			     nb_mesh2D_type  type);
void nb_mesh2D_copy(nb_mesh2D_t *mesh, const nb_mesh2D_t* srcmesh);
void nb_mesh2D_finish(nb_mesh2D_t *mesh);
nb_mesh2D_t* nb_mesh2D_create(nb_mesh2D_type type);
nb_mesh2D_t* nb_mesh2D_clone(nb_mesh2D_t* mesh);
void nb_mesh2D_clear(nb_mesh2D_t* mesh);
void nb_mesh2D_destroy(nb_mesh2D_t* mesh);

nb_mesh2D_type nb_mesh2D_get_type(const nb_mesh2D_t *mesh);

uint32_t nb_mesh2D_get_N_invtx(const nb_mesh2D_t *mesh);
uint32_t nb_mesh2D_get_N_insgm(const nb_mesh2D_t *mesh);
uint32_t nb_mesh2D_get_N_nodes(const nb_mesh2D_t *mesh);
uint32_t nb_mesh2D_get_N_edges(const nb_mesh2D_t *mesh);
uint32_t nb_mesh2D_get_N_elems(const nb_mesh2D_t *mesh);

double nb_mesh2D_node_get_x(const nb_mesh2D_t *mesh, uint32_t id);
double nb_mesh2D_node_get_y(const nb_mesh2D_t *mesh, uint32_t id);

uint32_t nb_mesh2D_edge_get_1n(const nb_mesh2D_t *mesh, uint32_t id);
uint32_t nb_mesh2D_edge_get_2n(const nb_mesh2D_t *mesh, uint32_t id);
void nb_mesh2D_edge_get_midpoint(const nb_mesh2D_t *mesh,
				 uint32_t face_id, double w,
				 double midpoint[2]);
double nb_mesh2D_edge_get_normal(const nb_mesh2D_t *mesh,
				 uint32_t face_id, double normal[2]);
double nb_mesh2D_edge_get_length(const nb_mesh2D_t *mesh,
				 uint32_t face_id);
double nb_mesh2D_elem_get_x(const nb_mesh2D_t *mesh, uint32_t id);
double nb_mesh2D_elem_get_y(const nb_mesh2D_t *mesh, uint32_t id);
double nb_mesh2D_elem_get_area(const nb_mesh2D_t *mesh, uint32_t id);
double nb_mesh2D_elem_get_radius(const nb_mesh2D_t *mesh, uint32_t id);
double nb_mesh2D_elem_get_apotem(const nb_mesh2D_t *mesh, uint32_t id);
uint32_t nb_mesh2D_elem_find_edge(const nb_mesh2D_t *mesh, uint32_t id,
				  uint16_t local_face_id);
double nb_mesh2D_elem_face_get_length(const nb_mesh2D_t *mesh,
				      uint32_t elem_id, uint16_t face_id);
double nb_mesh2D_elem_face_get_normal(const nb_mesh2D_t *mesh,
				      uint32_t elem_id, uint16_t face_id,
				      double normal[2]);
uint32_t nb_mesh2D_elem_face_get_left_ngb(const nb_mesh2D_t *mesh,
					  uint32_t elem_id,
					  uint16_t face_id);
uint32_t nb_mesh2D_elem_face_get_right_ngb(const nb_mesh2D_t *mesh,
					   uint32_t elem_id,
					   uint16_t face_id);
double nb_mesh2D_elem_ngb_get_normal(const nb_mesh2D_t *mesh,
				     uint32_t elem_id, uint16_t face_id,
				     double normal[2]);
uint16_t nb_mesh2D_elem_ngb_get_face(const nb_mesh2D_t *mesh,
				     uint32_t elem_id, uint32_t ngb_id);
uint32_t nb_mesh2D_elem_get_N_adj(const nb_mesh2D_t *mesh, uint32_t id);
uint32_t nb_mesh2D_elem_get_adj(const nb_mesh2D_t *mesh,
				uint32_t elem_id, uint8_t adj_id);
uint32_t nb_mesh2D_elem_get_ngb(const nb_mesh2D_t *mesh,
				uint32_t elem_id, uint8_t ngb_id);
bool nb_mesh2D_elem_has_ngb(const nb_mesh2D_t *mesh, uint32_t elem_id,
			    uint16_t face_id);
bool nb_mesh2D_elem_is_boundary(const nb_mesh2D_t *mesh, uint32_t elem_id);
uint32_t nb_mesh2D_get_invtx(const nb_mesh2D_t *mesh, uint32_t id);
uint32_t nb_mesh2D_insgm_get_N_nodes(const nb_mesh2D_t *mesh, uint32_t id);
uint32_t nb_mesh2D_insgm_get_N_subsgm(const nb_mesh2D_t *mesh, uint32_t id);
uint32_t nb_mesh2D_insgm_get_node(const nb_mesh2D_t *mesh,
				  uint32_t sgm_id, uint32_t node_id);
double nb_mesh2D_insgm_get_length(const nb_mesh2D_t *mesh, uint32_t sgm_id);
double nb_mesh2D_insgm_subsgm_get_length(const nb_mesh2D_t *mesh,
					 uint32_t sgm_id,
					 uint32_t subsgm_id);
void nb_mesh2D_insgm_get_elem_adj(const nb_mesh2D_t *mesh,
				  uint32_t **elem_adj);
void nb_mesh2D_load_graph(const nb_mesh2D_t *mesh, nb_graph_t *graph,
			  nb_mesh2D_graph_type type);

void nb_mesh2D_load_from_tessellator2D(nb_mesh2D_t *mesh, nb_tessellator2D_t *t2D);
void nb_mesh2D_set_nodal_permutation(nb_mesh2D_t *mesh,
				     const uint32_t *perm);
void nb_mesh2D_get_enveloping_box(const nb_mesh2D_t *mesh, double box[4]);
bool nb_mesh2D_is_vtx_inside(const nb_mesh2D_t *mesh, double x, double y);
void nb_mesh2D_extrapolate_elems_to_nodes(const nb_mesh2D_t *mesh,
					  uint8_t N_comp,
					  const double *elem_values,
					  double *nodal_values);
double nb_mesh2D_distort_with_field(nb_mesh2D_t *mesh,
				    nb_mesh2D_entity field_type,
				    double *disp, double max_disp);

void nb_mesh2D_build_model(const nb_mesh2D_t *mesh, nb_model_t *model);

/**
 * @brief Build a geometry model from a triangulation with disabled
 * triangles.
 * @param[in] mesh Input triangulation.
 * @param[in] elems_enabled Array of booleans indicating wich triangles are
 * enabled. NULL to enable all.
 * @param[out] model Model generated.
 * @param[out] N_input_vtx Stores the number of vertices forming
 * the boundary that appears in the original input.
 * @param[out] input_vtx Stores an array of the vertices ids in
 * the new model corresponding to original input vertices.
 */
void nb_mesh2D_build_model_disabled_elems(const nb_mesh2D_t *mesh,
					  const bool *elems_enabled,
					  nb_model_t *model,
					  uint32_t *N_input_vtx,
					  uint32_t **input_vtx);

#endif
