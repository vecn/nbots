#ifndef __NB_GEOMETRIC_BOT_MESH_PARTITION_H__
#define __NB_GEOMETRIC_BOT_MESH_PARTITION_H__

#include <stdint.h>

#include "nb/graph_bot.h"

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

double nb_partition_get_x_node(const nb_partition_t *part, uint32_t id);
double nb_partition_get_y_node(const nb_partition_t *part, uint32_t id);

uint32_t nb_partition_get_1n_edge(const nb_partition_t *part, uint32_t id);
uint32_t nb_partition_get_2n_edge(const nb_partition_t *part, uint32_t id);

double nb_partition_get_x_elem(const nb_partition_t *part, uint32_t id);
double nb_partition_get_y_elem(const nb_partition_t *part, uint32_t id);

uint32_t nb_partition_elem_get_N_adj(const nb_partition_t *part, uint32_t id);
uint32_t nb_partition_elem_get_adj(const nb_partition_t *part,
				   uint32_t elem_id, uint8_t adj_id);
uint32_t nb_partition_elem_get_N_ngb(const nb_partition_t *part, uint32_t id);
uint32_t nb_partition_elem_get_ngb(const nb_partition_t *part,
				   uint32_t elem_id, uint8_t ngb_id);

uint32_t nb_partition_get_invtx(const nb_partition_t *part, uint32_t id);
uint32_t nb_partition_get_N_nodes_x_insgm(const nb_partition_t *part,
					  uint32_t id);
uint32_t nb_partition_get_node_x_insgm(const nb_partition_t *part,
				       uint32_t sgm_id, uint32_t node_id);

void nb_partition_load_elem_graph(const nb_partition_t *part,
				  vcn_graph_t *graph);
void nb_partition_load_from_mesh(nb_partition_t *part,
				 const nb_mesh_t *const mesh);
void nb_partition_get_enveloping_box(const nb_partition_t *part,
				     double box[4]);

#endif
