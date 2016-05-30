#ifndef __NB_PDE_BOT_FINITE_ELEMENT_ELEMENT_STRUCT_H__
#define __NB_PDE_BOT_FINITE_ELEMENT_ELEMENT_STRUCT_H__

#include <stdint.h>

#include "nb/pde_bot/finite_element/element.h"

struct vcn_fem_elem_s {
	vcn_elem_id type;
	uint32_t N_nodes;
	uint32_t N_Gauss_points;
	double *psi, *eta;          /* Normalized space coordinates of Gauss points */
	double * gp_weight;         /* Integration weights of the Gauss points */
	double (**Ni)(double psi, double eta);       /* Shape functions */
	double (**dNi_dpsi)(double psi, double eta); /* Shape functions derivatives */
	double (**dNi_deta)(double psi, double eta); /* Shape functions derivatives */
	uint32_t (*get_closest_GP_to_ith_node)(uint32_t i);
};

#endif
