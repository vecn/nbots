#ifndef __NB_PDE_BOT_FINITE_ELEMENT_ELEMENT_H__
#define __NB_PDE_BOT_FINITE_ELEMENT_ELEMENT_H__

#include <stdint.h>

typedef struct vcn_fem_elem_s vcn_fem_elem_t;

typedef enum {
	NB_TRG_LINEAR,
} vcn_elem_id;

vcn_fem_elem_t* vcn_fem_elem_create(vcn_elem_id type);
void vcn_fem_elem_destroy(vcn_fem_elem_t* elemtype);
uint8_t vcn_fem_elem_get_N_gpoints(const vcn_fem_elem_t *const elemtype);
uint8_t vcn_fem_elem_get_N_nodes(const vcn_fem_elem_t *const elemtype);
double vcn_fem_elem_weight_gp(const vcn_fem_elem_t *const elemtype,
			      uint8_t gp_id);
double vcn_fem_elem_Ni(const vcn_fem_elem_t *const elemtype,
		       uint8_t node_id, uint8_t gp_id);
double vcn_fem_elem_dNi_dpsi(const vcn_fem_elem_t *const elemtype,
			     uint8_t node_id, uint8_t gp_id);
double vcn_fem_elem_dNi_deta(const vcn_fem_elem_t *const elemtype,
			     uint8_t node_id, uint8_t gp_id);

#endif
