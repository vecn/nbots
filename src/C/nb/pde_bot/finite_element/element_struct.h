#ifndef __NB_PDE_BOT_FINITE_ELEMENT_ELEMENT_STRUCT_H__
#define __NB_PDE_BOT_FINITE_ELEMENT_ELEMENT_STRUCT_H__

#include <stdint.h>

struct vcn_fem_elem_s {
	uint32_t N_nodes;
	uint32_t N_Gauss_points;
	double *psi, *eta;          /* Normalized space coordinates of Gauss points */
	double * gp_weight;         /* Integration weights of the Gauss points */
	double (**Ni)(double psi, double eta);       /* Shape functions */
	double (**dNi_dpsi)(double psi, double eta); /* Shape functions derivatives */
	double (**dNi_deta)(double psi, double eta); /* Shape functions derivatives */
};

#endif
