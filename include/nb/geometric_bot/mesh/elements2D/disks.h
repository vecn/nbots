#ifndef __NB_GEOMETRIC_BOT_MESH_ELEMENTS2D_DISKS_H__
#define __NB_GEOMETRIC_BOT_MESH_ELEMENTS2D_DISKS_H__

#include <stdint.h>
#include <stdbool.h>

#include "nb/geometric_bot/mesh/mesh2D.h"
#include "nb/graph_bot.h"

typedef struct nb_mshpack_s nb_mshpack_t;

uint32_t nb_mshpack_get_memsize(void);
void nb_mshpack_init(void *msh);
void nb_mshpack_finish(void *msh);
void nb_mshpack_copy(void *msh, const void *mshsrc);
void nb_mshpack_clear(void *msh);
uint32_t nb_mshpack_get_N_invtx(const void *msh);
uint32_t nb_mshpack_get_N_insgm(const void *msh);
uint32_t nb_mshpack_get_N_nodes(const void *msh);
uint32_t nb_mshpack_get_N_edges(const void *msh);
uint32_t nb_mshpack_get_N_elems(const void *msh);
double nb_mshpack_get_x_node(const void *msh, uint32_t id);
double nb_mshpack_get_y_node(const void *msh, uint32_t id);
uint32_t nb_mshpack_get_1n_edge(const void *msh, uint32_t id);
uint32_t nb_mshpack_get_2n_edge(const void *msh, uint32_t id);
double nb_mshpack_get_x_elem(const void *msh, uint32_t id);
double nb_mshpack_get_y_elem(const void *msh, uint32_t id);
double nb_mshpack_get_elem_radii(const void *msh, uint32_t id);
uint32_t nb_mshpack_elem_get_N_adj(const void *msh, uint32_t id);
uint32_t nb_mshpack_elem_get_adj(const void *msh,
				 uint32_t elem_id, uint8_t adj_id);
uint32_t nb_mshpack_elem_get_N_ngb(const void *msh, uint32_t id);
uint32_t nb_mshpack_elem_get_ngb(const void *msh,
				 uint32_t elem_id, uint8_t ngb_id);
uint32_t nb_mshpack_get_invtx(const void *msh, uint32_t id);
uint32_t nb_mshpack_get_N_nodes_x_insgm(const void *msh, uint32_t id);
uint32_t nb_mshpack_get_node_x_insgm(const void *msh, uint32_t sgm_id,
				     uint32_t node_id);
void nb_mshpack_load_elem_graph(const void *msh, nb_graph_t *graph);
void nb_mshpack_load_nodal_graph(const void *msh, nb_graph_t *graph);
void nb_mshpack_load_interelem_graph(const void *msh, nb_graph_t *graph);
void nb_mshpack_load_from_mesh_with_overlap(void *msh,
					    const nb_mesh_t *mesh,
					    double ov_factor);
void nb_mshpack_load_from_mesh(void *msh, const nb_mesh_t *mesh);
void nb_mshpack_get_enveloping_box(const void *msh, double box[4]);
bool nb_mshpack_is_vtx_inside(const void *msh, double x, double y);
void nb_mshpack_draw(const void *msh, const char *filename,
		     int width, int height);
void nb_mshpack_build_model(const void *msh, nb_model_t *model);
void nb_mshpack_build_model_disabled_elems(const void *msh,
					   const bool *elems_enabled,
					   nb_model_t *model,
					   uint32_t *N_input_vtx,
					   uint32_t **input_vtx);

#endif
