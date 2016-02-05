#ifndef __VCN_PDE_BOT_FINITE_ELEMENT_ELEMENT_H__
#define __VCN_PDE_BOT_FINITE_ELEMENT_ELEMENT_H__

typedef struct vcn_fem_elem_s vcn_fem_elem_t;

typedef enum {
	VCN_TRG_LINEAR,
} vcn_elem_id;

void vcn_fem_elem_create(vcn_elem_id id);
void vcn_fem_elem_destroy(vcn_fem_elem_t* elemtype);
uint32_t vcn_fem_elem_get_N_nodes(const vcn_fem_elem_t *const elemtype);
void* vcn_fem_elem_get_ith_shape_function
(const vcn_fem_elem_t *const elemtype, uint32_t i);
uint32_t vcn_fem_elem_get_closest_Gauss_Point_to_the_ith_node
(const vcn_fem_elem_t *const elemtype, uint32_t i);

#endif