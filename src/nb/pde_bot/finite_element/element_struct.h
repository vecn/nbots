#ifndef __NB_PDE_BOT_FINITE_ELEMENT_ELEMENT_STRUCT_H__
#define __NB_PDE_BOT_FINITE_ELEMENT_ELEMENT_STRUCT_H__

#include <stdint.h>

#include "nb/pde_bot/finite_element/element.h"

struct nb_fem_elem_s {
	nb_elem_id type;
	uint8_t N_nodes;
	uint8_t N_gp;
	double *gp_weight;    /* Integration weights of the Gauss pnt */
	double *Ni;           /* Shape functions at Gauss pnt */
	double *dNi_dpsi;     /* Shape functions derivatives at Gauss pnt */
	double *dNi_deta;     /* Shape functions derivatives at Gauss pnt */
};

#endif
