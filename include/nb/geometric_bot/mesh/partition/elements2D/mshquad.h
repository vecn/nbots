#ifndef __NB_GEOMETRIC_BOT_MESH_PARTITION_ELEMENTS2D_MSHQUAD_H__
#define __NB_GEOMETRIC_BOT_MESH_PARTITION_ELEMENTS2D_MSHQUAD_H__

#include <stdbool.h>
#include <stdint.h>

#include "nb/graph_bot.h"
#include "nb/geometric_bot/model/model2D.h"
#include "nb/geometric_bot/mesh/mesh2D.h"
#include "nb/geometric_bot/mesh/partition/info.h"

typedef struct nb_mshquad_s nb_mshquad_t;

uint32_t nb_mshquad_get_memsize(void);
void nb_mshquad_init(void *mshquad_ptr);
void nb_mshquad_copy(void *dest, const void *const src);
void nb_mshquad_finish(void *mshquad_ptr);

void* nb_mshquad_create(void);
void* nb_mshquad_clone(const void *const mshquad_ptr);
void nb_mshquad_destroy(void *mshquad_ptr);  
void nb_mshquad_clear(void *mshquad_ptr);

uint32_t nb_mshquad_get_N_invtx(const void *msh);
uint32_t nb_mshquad_get_N_insgm(const void *msh);
uint32_t nb_mshquad_get_N_nodes(const void *msh);
uint32_t nb_mshquad_get_N_edges(const void *msh);
uint32_t nb_mshquad_get_N_elems(const void *msh);
double nb_mshquad_node_get_x(const void *msh, uint32_t id);
double nb_mshquad_node_get_y(const void *msh, uint32_t id);
uint32_t nb_mshquad_edge_get_1n(const void *msh, uint32_t id);
uint32_t nb_mshquad_edge_get_2n(const void *msh, uint32_t id);
double nb_mshquad_elem_get_x(const void *msh, uint32_t id);
double nb_mshquad_elem_get_y(const void *msh, uint32_t id);
double nb_mshquad_elem_get_area(const void *msh, uint32_t id);
double nb_mshquad_elem_face_get_length(const void *msh,
				       uint32_t elem_id,
				       uint16_t face_id);
uint32_t nb_mshquad_elem_get_N_adj(const void *msh, uint32_t id);
uint32_t nb_mshquad_elem_get_adj(const void *msh,
				 uint32_t elem_id, uint8_t adj_id);
uint32_t nb_mshquad_elem_get_N_ngb(const void *msh, uint32_t id);
uint32_t nb_mshquad_elem_get_ngb(const void *msh,
				 uint32_t elem_id, uint8_t ngb_id);
bool nb_mshquad_elem_has_ngb(const void *msh, uint32_t elem_id,
			     uint16_t ngb_id);
bool nb_mshquad_elem_is_quad(const void *msh, uint32_t elem_id);
uint32_t nb_mshquad_get_invtx(const void *msh, uint32_t id);
uint32_t nb_mshquad_insgm_get_N_nodes(const void *msh, uint32_t id);
uint32_t nb_mshquad_insgm_get_node(const void *msh, uint32_t sgm_id,
				     uint32_t node_id);

void nb_mshquad_load_elem_graph(const void *mshquad,
				nb_graph_t *graph);
void nb_mshquad_load_nodal_graph(const void *mshquad,
				 nb_graph_t *graph);
void nb_mshquad_load_interelem_graph(const void *mshquad,
				     nb_graph_t *graph);

void nb_mshquad_load_from_mesh(void *mshquad, nb_mesh_t *mesh);

void nb_mshquad_get_enveloping_box(const void *msh, double box[4]);
bool nb_mshquad_is_vtx_inside(const void *msh, double x, double y);

double nb_mshquad_distort_with_field(void *msh, 
				     nb_partition_entity field_entity,
				     double *disp,
				     double max_disp);
void nb_mshquad_extrapolate_elems_to_nodes(const void *msh, uint8_t N_comp,
					   const double *elem_values,
					   double *nodal_values);
void nb_mshquad_build_model(const void *msh, nb_model_t *model);
void nb_mshquad_build_model_disabled_elems(const void *msh,
					   const bool *elems_enabled,
					   nb_model_t *model,
					   uint32_t *N_input_vtx,
					   uint32_t **input_vtx);

#endif