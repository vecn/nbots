#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "nb/math_bot.h"
#include "nb/cfreader_cat.h"
#include "nb/memory_bot.h"
#include "nb/container_bot.h"
#include "nb/solver_bot.h"
#include "nb/geometric_bot.h"
#include "nb/graph_bot.h"
#include "nb/pde_bot/material.h"
#include "nb/pde_bot/finite_element/element.h"
#include "nb/pde_bot/finite_element/gaussp_to_nodes.h"
#include "nb/pde_bot/common_solid_mechanics/analysis2D.h"
#include "nb/pde_bot/finite_element/solid_mechanics/static_elasticity2D.h"

#include "../utils.h"
#include "pipeline.h"
#include "nb/pde_bot/finite_element/solid_mechanics/static_plasticity2D.h"

#define POW2(a) ((a)*(a))

int compute_internal_forces
        (double *FI,
         double *stress,
         uint32_t N_elem,
         uint32_t N_nod,
         const nb_mesh2D_t *part,
         const nb_material_t *material,
         nb_analysis2D_t analysis2D,
         nb_analysis2D_params *params2D,
         const nb_fem_elem_t *const elem)
{
    int status = 1;
    N_nod = nb_fem_elem_get_N_nodes(elem);

	for (uint32_t i = 0; i < N_elem; i++) {

        int status_element = assemble_internal_forces_element(elem, i, stress, part, material, analysis2D, params2D, N_nod, FI);
		if (0 != status_element)
			goto EXIT;
	}
	status = 0;
	/* TEMPORAL */
	for (int k = 0; k < nb_mesh2D_get_N_nodes(part); k++){
        printf("FIx[%d]: %lf\t", k, FI[2*k]);
        printf("FIy[%d]: %lf\n", k, FI[2*k+1]);
	}
EXIT:
	return status;
}

int assemble_internal_forces_element(const nb_fem_elem_t *elem,
                uint32_t id, double *stress,
			    const nb_mesh2D_t *part,
			    const nb_material_t *material,
			    nb_analysis2D_t analysis2D,
			    nb_analysis2D_params *params2D,
			    uint32_t N_nod,
			    double *FI)
{
    double *FIe = nb_allocate_zero_mem(2*N_nod*sizeof(FIe));

    int status = integrate_elemental_FIe_vector
                    (elem, id, part, stress, params2D, N_nod, FIe);

	if (0 != status)
		goto CLEANUP;

    add_internal_forces_to_global_system
                    (elem, id, part, FI, FIe);
CLEANUP:
	nb_free_mem(FIe);

	return status;
}

int integrate_elemental_FIe_vector
             (const nb_fem_elem_t *elem, uint32_t id,
			  const nb_mesh2D_t *part, double *stress,
			  nb_analysis2D_params *params2D,
			  uint32_t N_nod,
			  double *FIe)
{
	int status = 1;
	/* Allocate Cartesian derivatives for each Gauss Point */
	uint32_t N_elem_nod = nb_fem_elem_get_N_nodes(elem);
	uint32_t deriv_memsize = N_elem_nod * sizeof(double);
	/*char *deriv_memblock = nb_allocate_mem(2 * deriv_memsize);
	double *dNi_dx = (void*) (deriv_memblock);
	double *dNi_dy = (void*) (deriv_memblock + deriv_memsize);*/
	double *dNi_dx = nb_soft_allocate_mem(deriv_memsize);
	memset(dNi_dx, 0, deriv_memsize);
	double *dNi_dy = nb_soft_allocate_mem(deriv_memsize);
    memset(dNi_dy, 0, deriv_memsize);

	uint8_t N_gp = nb_fem_elem_get_N_gpoints(elem);
	for (uint32_t j = 0; j < N_gp; j++) {
		double Jinv[4];
		double detJ = nb_fem_get_jacobian(elem, id, part, j, Jinv);

		if (nb_fem_elem_is_distorted(detJ))
			goto EXIT;

		nb_fem_get_derivatives(elem, j, Jinv, dNi_dx, dNi_dy);

        double thickness = params2D->thickness;
        internal_forces_sum_gauss_point(elem, id, j, thickness,
                                       detJ, part, dNi_dx, dNi_dy,
                                       stress, FIe);
	}
	status = 0;

EXIT:
	//nb_free_mem(2 * deriv_memsize, deriv_memblock);
	nb_soft_free_mem(deriv_memsize, dNi_dx);
	nb_soft_free_mem(deriv_memsize, dNi_dy);

	return status;
}

void internal_forces_sum_gauss_point(const nb_fem_elem_t *elem, uint32_t gp_id, uint8_t GP,
			      double thickness, double detJ, const nb_mesh2D_t *part,
			      double *dNi_dx, double *dNi_dy,
			      double *stress, double *FIe)
{
	/* Compute elemental stiffness matrix
	 *       _            _         _   _                   _                           _
	 *      |dNi/dx        |       | Sx  |                 |  (dNi/dx)*sx + (dNi/dy)*sxy |
	 * Bi = |       dNi/dy |   S = | Sy  |  Fi[i] = B'*S = |                             |
	 *      |dNi/dy dNi/dx_|       |_Sxy_|                 |_ (dNi/dy)*sy + (dNi/dx)*sx _|
	 */
	double wp = nb_fem_elem_weight_gp(elem, GP);

    uint8_t N_elem_node = nb_fem_elem_get_N_nodes(elem);
    memset(FIe, 0, 2*N_elem_node*sizeof(FIe));
	for (uint32_t i = 0; i < N_elem_node; i++) {
        uint32_t v1 = nb_mesh2D_elem_get_adj(part, gp_id, i);

		FIe[2*i] = (dNi_dx[i]*stress[3*gp_id] + dNi_dy[i]*stress[3*gp_id + 2]) *
		detJ * thickness * wp;
		FIe[2*i + 1] = (dNi_dy[i]*stress[3*gp_id + 1] + dNi_dx[i]*stress[3*gp_id]) *
		detJ * thickness * wp;
	}
}

void add_internal_forces_to_global_system(const nb_fem_elem_t *elem, uint32_t id,
				   const nb_mesh2D_t *part, double *FI,
				   double *FIe)
{
    uint8_t N_elem_node = nb_fem_elem_get_N_nodes(elem);
	for (uint32_t i = 0; i < N_elem_node; i++) {
		uint32_t v1 = nb_mesh2D_elem_get_adj(part, id, i);
		/* Add to global internal forces vector */
		FI[2*v1] += FIe[2*i];
		FI[2*v1 + 1] += FIe[2*i + 1];
	}
}
