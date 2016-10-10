#ifndef __NB_PDE_BOT_FINITE_ELEMENT_ELEMENT_H__
#define __NB_PDE_BOT_FINITE_ELEMENT_ELEMENT_H__

#include <stdint.h>

typedef struct nb_fem_elem_s nb_fem_elem_t;

typedef enum {
	NB_TRG_LINEAR,
	NB_QUAD_LINEAR
} nb_elem_id;

nb_fem_elem_t* nb_fem_elem_create(nb_elem_id type);
void nb_fem_elem_destroy(nb_fem_elem_t* elemtype);
uint8_t nb_fem_elem_get_N_gpoints(const nb_fem_elem_t *const elemtype);
uint8_t nb_fem_elem_get_N_nodes(const nb_fem_elem_t *const elemtype);
double nb_fem_elem_weight_gp(const nb_fem_elem_t *const elemtype,
			      uint8_t gp_id);
double nb_fem_elem_Ni(const nb_fem_elem_t *const elemtype,
		       uint8_t node_id, uint8_t gp_id);
double nb_fem_elem_dNi_dpsi(const nb_fem_elem_t *const elemtype,
			     uint8_t node_id, uint8_t gp_id);
double nb_fem_elem_dNi_deta(const nb_fem_elem_t *const elemtype,
			     uint8_t node_id, uint8_t gp_id);

#endif
