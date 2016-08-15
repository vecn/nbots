#include <stdint.h>

#include "nb/geometric_bot.h"
#include "nb/pde_bot/material.h"
#include "nb/pde_bot/finite_element/element.h"

#include "utils.h"

double nb_fem_get_jacobian(const vcn_fem_elem_t *elem, uint32_t id,
			     const nb_partition_t *part, int gp_id,
			     double Jinv[4])
{
	/* Compute Jacobian derivatives */
	double dx_dpsi = 0.0;
	double dy_dpsi = 0.0;
	double dx_deta = 0.0;
	double dy_deta = 0.0;

	uint8_t N_nodes = vcn_fem_elem_get_N_nodes(elem);
	for (uint32_t i = 0; i < N_nodes; i++) {
		uint32_t inode = nb_partition_elem_get_adj(part, id, i);
		double xi = nb_partition_get_x_node(part, inode);
		double yi = nb_partition_get_y_node(part, inode);
		double dNi_dpsi = vcn_fem_elem_dNi_dpsi(elem, i, gp_id);
		double dNi_deta = vcn_fem_elem_dNi_deta(elem, i, gp_id);
		dx_dpsi += dNi_dpsi * xi;
		dx_deta += dNi_deta * xi;
		dy_dpsi += dNi_dpsi * yi;
		dy_deta += dNi_deta * yi;
	}

	/* Compute Jacobian inverse and determinant */
	double detJ = dx_dpsi * dy_deta - dy_dpsi * dx_deta;
      
	Jinv[0] =  dy_deta / detJ;
	Jinv[1] = -dy_dpsi / detJ;
	Jinv[2] = -dx_deta / detJ;
	Jinv[3] =  dx_dpsi / detJ;

	return detJ;
}


bool nb_fem_elem_is_distorted(double detJ)
{
	return detJ < 0;
}

void nb_fem_get_derivatives(const vcn_fem_elem_t *elem,
			      int gp_id, double Jinv[4],
			      double *dNi_dx, double *dNi_dy)
{
	uint8_t N_nodes = vcn_fem_elem_get_N_nodes(elem);
	for (uint8_t i = 0; i < N_nodes; i++) {
		double dNi_dpsi = vcn_fem_elem_dNi_dpsi(elem, i, gp_id);
		double dNi_deta = vcn_fem_elem_dNi_deta(elem, i, gp_id);
		dNi_dx[i] = Jinv[0] * dNi_dpsi + Jinv[1] * dNi_deta;
		dNi_dy[i] = Jinv[2] * dNi_dpsi + Jinv[3] * dNi_deta;
	}
}
