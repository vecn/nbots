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
#include "nb/pde_bot/finite_element/solid_mechanics/static_elasticity2D.h"
#include "nb/pde_bot/finite_element/solid_mechanics/static_plasticity2D.h"
#include "nb/pde_bot/finite_element/solid_mechanics/plastic_stiffness_matrix.h"


#include "../utils.h"
#include "pipeline.h"

#define POW2(a) ((a)*(a))

int updated_stiffness_matrix(nb_sparse_t* K, nb_plastified_analysis2D *elem_regime,
                        const bool *elements_enabled, const nb_mesh2D_t *part,
                        const nb_fem_elem_t *elem,
                        const nb_material_t *material,
                        nb_analysis2D_t analysis2D,
                        nb_analysis2D_params *params2D,
                        double *elastic_strain,
                        uint32_t *simultaneous_elements,
                        double *aux_strain,
                        uint32_t *N_simultaneous_plastic_elem,
                        uint32_t *plastified_elem,
                        uint32_t *N_plastic_elem)
{
    int update_status;
    uint32_t sim_element;
    if(N_simultaneous_plastic_elem > 1) {
        for(int j = 0; j < N_simultaneous_plastic_elem; j++) {
            sim_element = simultaneous_elements[j];
            update_status = updated_plastified_stiffness_matrix(K, elem_regime[sim_element],
                                                                pipeline_elem_is_enabled(elements_enabled, sim_element),
                                                                part, sim_element, elem, material,
                                                                analysis2D, params2D);
            if(update_status != 0) {
                update_status = 4;
                return update_status;
            }
            elastic_strain[3*sim_element] = aux_strain[3*sim_element];
            elastic_strain[3*sim_element + 1] = aux_strain[3*sim_element+ 1];
            elastic_strain[3*sim_element + 2] = aux_strain[3*sim_element+ 2];
            elem_regime[sim_element] = NB_PLASTIC;
        }
    }
    else {
        update_status = updated_plastified_stiffness_matrix(K, elem_regime[plastified_elem[0]],
                                                            pipeline_elem_is_enabled(elements_enabled, sim_element), part,
                                                            plastified_elem[0], elem, material,
                                                            analysis2D, params2D);
        if(update_status != 0) {
            update_status = 4;
            return update_status;
        }

        elastic_strain[3*plastified_elem[0]] = aux_strain[3*plastified_elem[0]];
        elastic_strain[3*plastified_elem[0] + 1] = aux_strain[3*plastified_elem[0] + 1];
        elastic_strain[3*plastified_elem[0] + 2] = aux_strain[3*plastified_elem[0] + 2];
        elem_regime[plastified_elem[0]] = NB_PLASTIC;
    }
    N_plastic_elem[0] -= 1;

    return update_status;
}

int updated_plastified_stiffness_matrix(nb_sparse_t* K, nb_plastified_analysis2D elem_regime,
                                    bool is_enabled, const nb_mesh2D_t *part,
                                    uint32_t plastified_elem,
                                    const nb_fem_elem_t *elem,
                                    const nb_material_t *material,
                                    nb_analysis2D_t analysis2D,
                                    nb_analysis2D_params *params2D)
 {
    int status = 1;
	uint32_t N_elem = nb_mesh2D_get_N_elems(part);

        /*printf("%s\n", is_enabled ? "true" : "false");*/
		int status_element =
			assemble_plastified_element(elem, plastified_elem, part, material, is_enabled,
					 elem_regime, analysis2D, params2D, K);

		if (0 != status_element)
			goto EXIT;

	status = 0;
EXIT:
	return status;
}

static int assemble_plastified_element(const nb_fem_elem_t *elem, uint32_t id,
			    const nb_mesh2D_t *part,
			    const nb_material_t *material,
			    bool is_enabled,
			    nb_plastified_analysis2D elem_regime,
			    nb_analysis2D_t analysis2D,
			    nb_analysis2D_params *params2D,
			    nb_sparse_t *K)
{
	double De[4] = {1e-6, 1e-6, 1e-6, 1e-6};
	double Dp[4] = {1e-6, 1e-6, 1e-6, 1e-6};
	double density = 1e-6;
	if (is_enabled) {
        nb_pde_get_plastified_constitutive_matrix(Dp, material, analysis2D, elem_regime);
        nb_pde_get_constitutive_matrix(De, material, analysis2D);
	}

	uint8_t N_nodes = nb_fem_elem_get_N_nodes(elem);
	double *Ke_elastic = nb_soft_allocate_mem(4 * POW2(N_nodes) * sizeof(*Ke_elastic));
	int status_elastic = integrate_plastic_elemental_system(elem, id, De, part, params2D, Ke_elastic);
	if (0 != status_elastic)
		goto CLEANUP;

    double *Ke_plastic = nb_soft_allocate_mem(4 * POW2(N_nodes) * sizeof(*Ke_elastic));
    int status_plastic = integrate_plastic_elemental_system(elem, id, Dp, part, params2D, Ke_plastic);
	if (0 != status_plastic)
		goto CLEANUP;

    modify_global_system(elem, id, part, Ke_elastic, Ke_plastic, K);

CLEANUP:
	nb_soft_free_mem(4 * POW2(N_nodes) * sizeof(*Ke_elastic), Ke_elastic);
	nb_soft_free_mem(4 * POW2(N_nodes) * sizeof(*Ke_plastic), Ke_plastic);

	return status_plastic;
}

static int integrate_plastic_elemental_system
            (const nb_fem_elem_t *elem, uint32_t id,
			 double D[4], const nb_mesh2D_t *part,
			 nb_analysis2D_params *params2D,
			 double *Ke)
{
	int status = 1;

	uint8_t N_nodes = nb_fem_elem_get_N_nodes(elem);
	memset(Ke, 0, 4 * POW2(N_nodes) * sizeof(*Ke));

	/* Allocate Cartesian derivatives for each Gauss Point */
	uint32_t deriv_memsize = N_nodes * sizeof(double);
	char *deriv_memblock = nb_soft_allocate_mem(2 * deriv_memsize);
	double *dNi_dx = (void*) (deriv_memblock);
	double *dNi_dy = (void*) (deriv_memblock + deriv_memsize);

	uint8_t N_gp = nb_fem_elem_get_N_gpoints(elem);
	for (uint32_t j = 0; j < N_gp; j++) {
		double Jinv[4];
		double detJ = nb_fem_get_jacobian(elem, id, part, j, Jinv);

		if (nb_fem_elem_is_distorted(detJ))
			goto EXIT;

		nb_fem_get_derivatives(elem, j, Jinv, dNi_dx, dNi_dy);

		double thickness = params2D->thickness;
		sum_plastic_gauss_point(elem, j, D, thickness,
					 detJ, dNi_dx, dNi_dy, Ke);
	}
	status = 0;
EXIT:
	nb_soft_free_mem(2 * deriv_memsize, deriv_memblock);
	return status;
}

void sum_plastic_gauss_point(const nb_fem_elem_t *elem, int gp_id,
			      double D[4], double thickness, double detJ,
			      double *dNi_dx, double *dNi_dy,
			      double *Ke)
{
	/* Compute elemental stiffness matrix
	 *       _            _         _        _        _  _
	 *      |dNi/dx        |       | D0 D1    |      |  |
	 * Bi = |       dNi/dy |   D = | D1 D2    |  K = |  | B'DB dx dy
	 *      |dNi/dy dNi/dx_|       |_      D3_|     _| _|

	 * B  = [B1 B2...Bi... Bn]
	 *
	 */
	double wp = nb_fem_elem_weight_gp(elem, gp_id);

	uint8_t N_nodes = nb_fem_elem_get_N_nodes(elem);
	for (uint32_t i = 0; i < N_nodes; i++) {
		double Ni = nb_fem_elem_Ni(elem, i, gp_id);
		for (uint32_t j = 0; j < N_nodes; j++) {
			/*  Integrating elemental siffness matrix */
			Ke[(i * 2)*(2 * N_nodes) + (j * 2)] +=
				(dNi_dx[i]*dNi_dx[j]*D[0] +
				 dNi_dy[i]*dNi_dy[j]*D[3]) *
				detJ * thickness * wp;
			Ke[(i * 2)*(2 * N_nodes) + (j*2+1)] +=
				(dNi_dx[i]*dNi_dy[j]*D[1] +
				 dNi_dy[i]*dNi_dx[j]*D[3]) *
				detJ * thickness * wp;
			Ke[(i*2+1)*(2 * N_nodes) + (j * 2)] +=
				(dNi_dy[i]*dNi_dx[j]*D[1] +
				 dNi_dx[i]*dNi_dy[j]*D[3]) *
				detJ * thickness * wp;
			Ke[(i*2+1)*(2 * N_nodes) + (j*2+1)] +=
				(dNi_dy[i]*dNi_dy[j]*D[2] +
				 dNi_dx[i]*dNi_dx[j]*D[3]) *
				detJ * thickness * wp;
        }
    }
}

void modify_global_system(const nb_fem_elem_t *elem, uint32_t id, const nb_mesh2D_t *part,
                          double *Ke_elastic, double *Ke_plastic, nb_sparse_t *K)
{
	uint8_t N_nodes = nb_fem_elem_get_N_nodes(elem);
	for (uint32_t i = 0; i < N_nodes; i++) {
		uint32_t v1 = nb_mesh2D_elem_get_adj(part, id, i);
		for (uint32_t j = 0; j < N_nodes; j++) {
			uint32_t v2 = nb_mesh2D_elem_get_adj(part, id, j);
			nb_sparse_add(K, v1 * 2, v2 * 2,
				       Ke_plastic[(i * 2)*(2 * N_nodes) +
					  (j * 2)]);
			nb_sparse_add(K, v1*2 + 1, v2 * 2,
				       Ke_plastic[(i*2 + 1)*(2 * N_nodes) +
					  (j * 2)]);
			nb_sparse_add(K, v1 * 2, v2*2 + 1,
				       Ke_plastic[(i * 2)*(2 * N_nodes) +
					  (j*2+1)]);
			nb_sparse_add(K, v1*2 + 1, v2*2 + 1,
				       Ke_plastic[(i*2 + 1)*(2 * N_nodes) +
					  (j*2+1)]);
			nb_sparse_substract(K, v1 * 2, v2 * 2,
				       Ke_elastic[(i * 2)*(2 * N_nodes) +
					  (j * 2)]);
			nb_sparse_substract(K, v1*2 + 1, v2 * 2,
				       Ke_elastic[(i*2 + 1)*(2 * N_nodes) +
					  (j * 2)]);
			nb_sparse_substract(K, v1 * 2, v2*2 + 1,
				       Ke_elastic[(i * 2)*(2 * N_nodes) +
					  (j*2+1)]);
			nb_sparse_substract(K, v1*2 + 1, v2*2 + 1,
				       Ke_elastic[(i*2 + 1)*(2 * N_nodes) +
					  (j*2+1)]);
		}
	}
}
