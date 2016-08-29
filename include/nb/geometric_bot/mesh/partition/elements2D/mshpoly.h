#ifndef __NB_GEOMETRIC_BOT_MESH_PARTITION_ELEMENTS2D_MSHPOLY_H__
#define __NB_GEOMETRIC_BOT_MESH_PARTITION_ELEMENTS2D_MSHPOLY_H__

#include <stdbool.h>
#include <stdint.h>

#include "nb/graph_bot.h"
#include "nb/geometric_bot/model/model2D.h"
#include "nb/geometric_bot/mesh/mesh2D.h"
#include "nb/geometric_bot/mesh/partition/info.h"

typedef struct nb_mshpoly_s nb_mshpoly_t;

uint32_t nb_mshpoly_get_memsize(void);
void nb_mshpoly_init(void *mshpoly_ptr);
void nb_mshpoly_copy(void *dest, const void *const src);
void nb_mshpoly_finish(void *mshpoly_ptr);

void* nb_mshpoly_create(void);
void* nb_mshpoly_clone(const void *const mshpoly_ptr);
void nb_mshpoly_destroy(void *mshpoly_ptr);  

void nb_mshpoly_clear(void *mshpoly_ptr);

uint32_t nb_mshpoly_get_N_invtx(const void *msh);
uint32_t nb_mshpoly_get_N_insgm(const void *msh);
uint32_t nb_mshpoly_get_N_nodes(const void *msh);
uint32_t nb_mshpoly_get_N_edges(const void *msh);
uint32_t nb_mshpoly_get_N_elems(const void *msh);
double nb_mshpoly_node_get_x(const void *msh, uint32_t id);
double nb_mshpoly_node_get_y(const void *msh, uint32_t id);
uint32_t nb_mshpoly_edge_get_1n(const void *msh, uint32_t id);
uint32_t nb_mshpoly_edge_get_2n(const void *msh, uint32_t id);
double nb_mshpoly_elem_get_x(const void *msh, uint32_t id);
double nb_mshpoly_elem_get_y(const void *msh, uint32_t id);
double nb_mshpoly_elem_get_area(const void *msh, uint32_t id);
double nb_mshpoly_elem_face_get_length(const void *msh,
				       uint32_t elem_id,
				       uint16_t face_id);
uint32_t nb_mshpoly_elem_get_N_adj(const void *msh, uint32_t id);
uint32_t nb_mshpoly_elem_get_adj(const void *msh,
				 uint32_t elem_id, uint8_t adj_id);
uint32_t nb_mshpoly_elem_get_N_ngb(const void *msh, uint32_t id);
uint32_t nb_mshpoly_elem_get_ngb(const void *msh,
				 uint32_t elem_id, uint8_t ngb_id);
bool nb_mshpoly_elem_has_ngb(const void *msh, uint32_t elem_id,
			     uint16_t ngb_id);
uint32_t nb_mshpoly_get_invtx(const void *msh, uint32_t id);
uint32_t nb_mshpoly_insgm_get_N_nodes(const void *msh, uint32_t id);
uint32_t nb_mshpoly_insgm_get_node(const void *msh, uint32_t sgm_id,
				     uint32_t node_id);

void nb_mshpoly_load_elem_graph(const void *mshpoly,
				nb_graph_t *graph);
void nb_mshpoly_load_nodal_graph(const void *mshpoly,
				 nb_graph_t *graph);
void nb_mshpoly_load_interelem_graph(const void *mshpoly,
				      nb_graph_t *graph);

void nb_mshpoly_get_enveloping_box(const void *msh, double box[4]);
bool nb_mshpoly_is_vtx_inside(const void *msh, double x, double y);

double nb_mshpoly_distort_with_field(void *msh, 
				     nb_partition_entity field_entity,
				     double *disp,
				     double max_disp);
void nb_mshpoly_extrapolate_elems_to_nodes(const void *msh, uint8_t N_comp,
					   const double *elem_values,
					   double *nodal_values);

void nb_mshpoly_build_model(const void *msh, nb_model_t *model);
void nb_mshpoly_build_model_disabled_elems(const void *msh,
					   const bool *elems_enabled,
					   nb_model_t *model,
					   uint32_t *N_input_vtx,
					   uint32_t **input_vtx);

void nb_mshpoly_load_from_mesh(void *mshpoly, nb_mesh_t *mesh);

void nb_mshpoly_Lloyd_iteration(void *mshpoly, uint32_t max_iter,
				double (*density)(const double[2],
						  const void *data),
				const void *density_data);

#endif
