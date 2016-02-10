#include <stdlib.h>

#include "nb/pde_bot/finite_element/element.h"
#include "nb/pde_bot/finite_element/elements/trg_linear.h"

#include "../element_struct.h"

#define INV_3 (0.33333333333333333333333333333)

static double N1(double psi, double eta);
static double N2(double psi, double eta);
static double N3(double psi, double eta);
static double dN1_dpsi(double psi, double eta);
static double dN1_deta(double psi, double eta);
static double dN2_dpsi(double psi, double eta);
static double dN2_deta(double psi, double eta);
static double dN3_dpsi(double psi, double eta);
static double dN3_deta(double psi, double eta);

void vcn_fem_elem_init_trg_linear(vcn_fem_elem_t *elem)
{
	elem->N_nodes = 3;
	elem->N_Gauss_points = 1;

	elem->Ni = calloc(elem->N_nodes, sizeof(*(elem->Ni)));
	elem->dNi_dpsi = calloc(elem->N_nodes,
				    sizeof(*(elem->dNi_dpsi)));
	elem->dNi_deta = calloc(elem->N_nodes,
				    sizeof(*(elem->dNi_deta)));  
	elem->Ni[0] = N1;
	elem->Ni[1] = N2;
	elem->Ni[2] = N3;
	elem->dNi_dpsi[0] = dN1_dpsi;
	elem->dNi_dpsi[1] = dN2_dpsi;
	elem->dNi_dpsi[2] = dN3_dpsi;
	elem->dNi_deta[0] = dN1_deta;
	elem->dNi_deta[1] = dN2_deta;
	elem->dNi_deta[2] = dN3_deta;
  
	elem->psi = malloc(elem->N_Gauss_points * sizeof(*(elem->psi)));
	elem->eta = malloc(elem->N_Gauss_points * sizeof(*(elem->eta)));
	elem->gp_weight = malloc(elem->N_Gauss_points *
				 sizeof(*(elem->gp_weight)));

	elem->psi[0] = INV_3;
	elem->eta[0] = INV_3;

	elem->gp_weight[0] = 0.5;
}

static inline double N1(double psi, double eta){
	return 1.0 - psi - eta;
}

static inline double N2(double psi, double eta){
	return psi;
}

static inline double N3(double psi, double eta){
	return eta;
}

static inline double dN1_dpsi(double psi, double eta){
	return -1.0;
}

static inline double dN1_deta(double psi, double eta){
	return -1.0;
}

static inline double dN2_dpsi(double psi, double eta){
	return 1.0;
}

static inline double dN2_deta(double psi, double eta){
	return 0.0;
}

static inline double dN3_dpsi(double psi, double eta){
	return 0.0;
}

static inline double dN3_deta(double psi, double eta){
	return 1.0;
}
