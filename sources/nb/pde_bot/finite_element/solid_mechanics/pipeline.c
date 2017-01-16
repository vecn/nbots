#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "nb/math_bot.h"
#include "nb/memory_bot.h"
#include "nb/container_bot.h"
#include "nb/solver_bot.h"
#include "nb/geometric_bot.h"
#include "nb/graph_bot.h"
#include "nb/pde_bot/material.h"
#include "nb/pde_bot/finite_element/element.h"
#include "nb/pde_bot/finite_element/gaussp_to_nodes.h"
#include "nb/pde_bot/finite_element/solid_mechanics/static_elasticity2D.h"

#include "../utils.h"
#include "pipeline.h"

#define POW2(a) ((a)*(a))

static int assemble_element(const nb_fem_elem_t *elem, uint32_t id,
			    const nb_mesh2D_t *part,
			    const nb_material_t *material,
			    bool is_enabled,
			    nb_analysis2D_t analysis2D,
			    nb_analysis2D_params *params2D,
			    bool enable_self_weight,
			    double gravity[2],
			    nb_sparse_t *K, double *M, double *F);
static int integrate_elemental_system
		       	(const nb_fem_elem_t *elem, uint32_t id,
			 double D[4], double density, double gravity[2],
			 const nb_mesh2D_t *part,
			 nb_analysis2D_params *params2D,
			 bool enable_self_weight,
			 double *Ke, double *Me, double *Fe);
static int get_element_strain(uint32_t id, double *strain,
			      const nb_mesh2D_t *const part,
			      double *displacement,
			      const nb_fem_elem_t *const elem);

int pipeline_assemble_system
		(nb_sparse_t* K, double* M, double *F,
		 const nb_mesh2D_t *const part,
		 const nb_fem_elem_t *const elem,
		 const nb_material_t *const material,
		 bool enable_self_weight,
		 double gravity[2],
		 nb_analysis2D_t analysis2D,
		 nb_analysis2D_params *params2D,
		 const bool* elements_enabled /* NULL to enable all */)
{
	int status = 1;
	uint32_t N_elem = nb_mesh2D_get_N_elems(part);
	nb_sparse_reset(K);
	if (NULL != M)
		memset(M, 0, nb_sparse_get_size(K) * sizeof(*M));
	memset(F, 0, nb_sparse_get_size(K) * sizeof(*F));

	for (uint32_t i = 0; i < N_elem; i++) {
		bool is_enabled = pipeline_elem_is_enabled(elements_enabled, i);
        /*printf("%s\n", is_enabled ? "true" : "false");*/
		int status_element =
			assemble_element(elem, i, part, material, is_enabled,
					 analysis2D, params2D,
					 enable_self_weight, gravity,
					 K, M, F);
		if (0 != status_element)
			goto EXIT;
	}
	status = 0;
EXIT:
	return status;
}

bool pipeline_elem_is_enabled(const bool *elements_enabled, uint32_t id)
{
	bool is_enabled = true;
	if (NULL != elements_enabled)
		is_enabled = elements_enabled[id];
	return is_enabled;
}

static int assemble_element(const nb_fem_elem_t *elem, uint32_t id,
			    const nb_mesh2D_t *part,
			    const nb_material_t *material,
			    bool is_enabled,
			    nb_analysis2D_t analysis2D,
			    nb_analysis2D_params *params2D,
			    bool enable_self_weight,
			    double gravity[2],
			    nb_sparse_t *K, double *M, double *F)
{
	double D[4] = {1e-6, 1e-6, 1e-6, 1e-6};
	double density = 1e-6;
	if (is_enabled) {
		nb_pde_get_constitutive_matrix(D, material, analysis2D);
		density = nb_material_get_density(material);
	}

	uint8_t N_nodes = nb_fem_elem_get_N_nodes(elem);
	double* Ke = nb_allocate_mem(4 * POW2(N_nodes) * sizeof(*Ke));
	double* Me = NULL;
	if(M != NULL)
		Me = nb_allocate_mem(2 * N_nodes * sizeof(*Me));
	double* Fe = nb_allocate_mem(2 * N_nodes * sizeof(*Fe));

	int status = integrate_elemental_system(elem, id, D, density,
						gravity, part, params2D,
						enable_self_weight,
						Ke, Me, Fe);
	if (0 != status)
		goto CLEANUP;

	pipeline_add_to_global_system(elem, id, part, Ke, Me, Fe, K, M, F);

CLEANUP:
	nb_free_mem(Ke);
	if (NULL != M)
		nb_free_mem(Me);
	nb_free_mem(Fe);
	return status;
}

static int integrate_elemental_system
            (const nb_fem_elem_t *elem, uint32_t id,
			 double D[4], double density, double gravity[2],
			 const nb_mesh2D_t *part,
			 nb_analysis2D_params *params2D,
			 bool enable_self_weight,
			 double *Ke, double *Me, double *Fe)
{
	int status = 1;

	double fx = 0.0;
	double fy = 0.0;
	if (enable_self_weight) {
		fx = gravity[0] * density;
		fy = gravity[1] * density;
	}

	uint8_t N_nodes = nb_fem_elem_get_N_nodes(elem);
	memset(Ke, 0, 4 * POW2(N_nodes) * sizeof(*Ke));
	if(Me != NULL)
		memset(Me, 0, 2 * N_nodes * sizeof(*Me));
	memset(Fe, 0, 2 * N_nodes * sizeof(*Fe));

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
		pipeline_sum_gauss_point(elem, j, D, density, thickness,
					 detJ, dNi_dx, dNi_dy, fx, fy,
					 Ke, Me, Fe);
	}
	status = 0;
EXIT:
	nb_soft_free_mem(2 * deriv_memsize, deriv_memblock);
	return status;
}

void pipeline_sum_gauss_point(const nb_fem_elem_t *elem, int gp_id,
			      double D[4], double density,
			      double thickness, double detJ,
			      double *dNi_dx, double *dNi_dy,
			      double fx, double fy,
			      double *Ke, double *Me, double *Fe)
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

			/*  Integrating elemental mass matrix */
			double Nj = nb_fem_elem_Ni(elem, j, gp_id);
			if (NULL != Me) {
				/* OPPORTUNITY: Allocate one for each comp */
				double integral = Ni * Nj * density * detJ *
					thickness * wp;
				Me[i * 2] += integral;
				Me[i*2+1] += integral;
			}
		}
		/* Compute elemental forces vector */
		/* OPPORTUNITY: Allocate just one for each component */
		double integral = Ni * detJ * thickness * wp;
		Fe[i * 2] += integral * fx;
		Fe[i*2+1] += integral * fy;
	}
}

void pipeline_add_to_global_system(const nb_fem_elem_t *elem, uint32_t id,
				   const nb_mesh2D_t *part,
				   double *Ke, double *Me, double *Fe,
				   nb_sparse_t *K, double *M, double *F)
{
	uint8_t N_nodes = nb_fem_elem_get_N_nodes(elem);
	for (uint32_t i = 0; i < N_nodes; i++) {
		uint32_t v1 = nb_mesh2D_elem_get_adj(part, id, i);
		for (uint32_t j = 0; j < N_nodes; j++) {
			uint32_t v2 = nb_mesh2D_elem_get_adj(part, id, j);
			nb_sparse_add(K, v1 * 2, v2 * 2,
				       Ke[(i * 2)*(2 * N_nodes) +
					  (j * 2)]);
			nb_sparse_add(K, v1*2 + 1, v2 * 2,
				       Ke[(i*2 + 1)*(2 * N_nodes) +
					  (j * 2)]);
			nb_sparse_add(K, v1 * 2, v2*2 + 1,
				       Ke[(i * 2)*(2 * N_nodes) +
					  (j*2+1)]);
			nb_sparse_add(K, v1*2 + 1, v2*2 + 1,
				       Ke[(i*2 + 1)*(2 * N_nodes) +
					  (j*2+1)]);
		}
		/* Add to global mass matrix */
		if (NULL != M) {
			M[v1 * 2] += Me[i * 2];
			M[v1*2+1] += Me[i*2+1];
		}
		/* Add to global forces vector */
		F[v1 * 2] += Fe[i * 2];
		F[v1*2+1] += Fe[i*2+1];
	}
}

void pipeline_compute_strain(double *strain,
			     const nb_mesh2D_t *const part,
			     double *displacement,
			     const nb_fem_elem_t *const elem)
{
	uint32_t N_elem = nb_mesh2D_get_N_elems(part);

	uint8_t N_gp = nb_fem_elem_get_N_gpoints(elem);
	memset(strain, 0, 3 * N_gp * N_elem * sizeof(*strain));
	for (uint32_t i = 0; i < N_elem; i++)
		get_element_strain(i, strain, part, displacement, elem);
}

static int get_element_strain(uint32_t id, double *strain,
			      const nb_mesh2D_t *const part,
			      double *displacement,
			      const nb_fem_elem_t *const elem)
{
	int status = 1;

	/* Allocate Cartesian derivatives for each Gauss Point */
	uint8_t N_nodes = nb_fem_elem_get_N_nodes(elem);
	uint32_t deriv_memsize = N_nodes * sizeof(double);
	char *deriv_memblock = nb_soft_allocate_mem(2 * deriv_memsize);
	double *dNi_dx = (void*) (deriv_memblock);
	double *dNi_dy = (void*) (deriv_memblock + deriv_memsize);

	//printf("N nodes: %d\n", N_nodes);
	//printf("Deriv memsize: %d\n", deriv_memsize);
	/* Integrate domain */
	uint8_t N_gp = nb_fem_elem_get_N_gpoints(elem);
	//printf("N gauss points: %d\n", N_gp);
	for (uint32_t j = 0; j < N_gp; j++) {
		double Jinv[4];
		double detJ = nb_fem_get_jacobian(elem, id, part, j, Jinv);
        //printf("Det J: %lf\n", detJ);

       // printf("Jinv[%d]: %lf\t Jinv[%d]: %lf\t Jinv[%d]: %lf\t Jinv[%d]: %lf\n", j, Jinv[0], j, Jinv[1], j, Jinv[2], j, Jinv[3]);

		if (nb_fem_elem_is_distorted(detJ)) {
            printf("Element %d is too distorted to be computed.\n", id);
            goto EXIT;
		}

		nb_fem_get_derivatives(elem, j, Jinv, dNi_dx, dNi_dy);

		for(int k = 0; k < N_nodes; k++){
           // printf("dNi_dx[%d, %d]: %lf\t", j, k, dNi_dx[k]);
           // printf("dNi_dy[%d, %d]: %lf\n", j, k, dNi_dy[k]);
		}

		uint32_t idx = id * N_gp + j;
		//printf("Idx: %d\n", idx);

		/* Compute Strain at Gauss Point */
		for (uint32_t i = 0; i < N_nodes; i++) {
			uint32_t inode = nb_mesh2D_elem_get_adj(part, id, i);
			//printf("I node: %d\n", inode);
			//printf("Dx[%d]: %lf\t", i, displacement[2*inode]);
			//printf("Dy[%d]: %lf\n", i, displacement[2*inode +1]);
			strain[idx * 3] += dNi_dx[i] * displacement[inode * 2];
			strain[idx*3+1] += dNi_dy[i] * displacement[inode*2+1];
			strain[idx*3+2] += (dNi_dy[i] * displacement[inode * 2] +
					    dNi_dx[i] * displacement[inode*2+1]);
           // printf("Ex: %lf\t", strain[3*idx]);
            //printf("Ey: %lf\t", strain[3*idx +1]);
            //printf("Exy: %lf\n", strain[3*idx +2]);
		}
	}
	status = 0;
EXIT:
	nb_soft_free_mem(2 * deriv_memsize, deriv_memblock);
	return status;
}

void pipeline_compute_main_stress(double *stress,
				  double *main_stress,
				  uint32_t N_elements,
				  const nb_fem_elem_t *const elem)
{
	uint8_t N_gp = nb_fem_elem_get_N_gpoints(elem);
	for(uint32_t i=0; i < N_elements; i++){
		for(uint32_t j=0; j < N_gp; j++){
			uint32_t idx = i * N_gp + j;
			double sigma_avg = 0.5 * (stress[idx * 3] +
						  stress[idx*3+1]);
			double R = sqrt(0.25 *
					POW2(stress[idx * 3] -
					     stress[idx*3+1]) +
					POW2(stress[idx*3+2]));
			main_stress[idx * 2] = sigma_avg + R;
			main_stress[idx*2+1] = sigma_avg - R;
		}
	}
}
