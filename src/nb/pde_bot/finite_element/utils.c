#include <stdint.h>

#include "nb/geometric_bot.h"
#include "nb/pde_bot/material.h"
#include "nb/pde_bot/finite_element/element.h"


#include "utils.h"
#include "element_struct.h"

double nb_fem_get_jacobian(const vcn_fem_elem_t *elem, uint32_t id,
			     const vcn_msh3trg_t *mesh, int gp_id,
			     double Jinv[4])
{
	uint32_t *conn_mtx = mesh->vertices_forming_triangles;

	/* Compute Jacobian derivatives */
	double dx_dpsi = 0.0;
	double dy_dpsi = 0.0;
	double dx_deta = 0.0;
	double dy_deta = 0.0;
	for (uint32_t i = 0; i < elem->N_nodes; i++) {
		uint32_t inode = conn_mtx[id * elem->N_nodes + i];
		double xi = mesh->vertices[inode * 2];
		double yi = mesh->vertices[inode*2+1];
		dx_dpsi += elem->dNi_dpsi[i](elem->psi[gp_id],
					     elem->eta[gp_id]) * xi;
		dx_deta += elem->dNi_deta[i](elem->psi[gp_id],
					     elem->eta[gp_id]) * xi;
		dy_dpsi += elem->dNi_dpsi[i](elem->psi[gp_id],
					     elem->eta[gp_id]) * yi;
		dy_deta += elem->dNi_deta[i](elem->psi[gp_id],
					     elem->eta[gp_id]) * yi;
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
	for (uint32_t i = 0; i < elem->N_nodes; i++) {
		dNi_dx[i] = Jinv[0]*elem->dNi_dpsi[i](elem->psi[gp_id],
						      elem->eta[gp_id]) + 
			Jinv[1]*elem->dNi_deta[i](elem->psi[gp_id],
						  elem->eta[gp_id]);
		dNi_dy[i] = Jinv[2]*elem->dNi_dpsi[i](elem->psi[gp_id],
						      elem->eta[gp_id]) + 
			Jinv[3]*elem->dNi_deta[i](elem->psi[gp_id],
						  elem->eta[gp_id]);
	}
}
