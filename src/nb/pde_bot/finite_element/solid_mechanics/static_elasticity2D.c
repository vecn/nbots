#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "nb/math_bot.h"
#include "nb/cfreader_cat.h"
#include "nb/eigen_bot.h"
#include "nb/container_bot.h"
#include "nb/graph_bot.h"
#include "nb/pde_bot/material.h"
#include "nb/pde_bot/boundary_conditions/bcond.h"
#include "nb/pde_bot/boundary_conditions/bcond_iter.h"
#include "nb/pde_bot/finite_element/element.h"
#include "nb/pde_bot/finite_element/gaussp_to_nodes.h"
#include "nb/pde_bot/finite_element/solid_mechanics/analysis2D.h"
#include "nb/pde_bot/finite_element/solid_mechanics/static_elasticity2D.h"
  
#include "../element_struct.h"
#include "pipeline.h"

#define POW2(a) ((a)*(a))

static int solver(const vcn_sparse_t *const A,
		  const double *const b, double* x);

int vcn_fem_compute_2D_Solid_Mechanics
			(const vcn_msh3trg_t *const mesh,
			 const vcn_fem_elem_t *const elemtype,
			 const vcn_fem_material_t *const material,
			 const nb_bcond_t *const bcond,
			 char enable_self_weight,
			 double gravity[2],
			 nb_analysis2D_t analysis2D,
			 nb_analysis2D_params *params2D,
			 const bool *elements_enabled, /* NULL to enable all */
			 double *displacement, /* Output */
			 double *strain       /* Output */)
{
	int status = 1;
	vcn_graph_t *graph = vcn_msh3trg_create_vtx_graph(mesh);
	vcn_sparse_t *K = vcn_sparse_create(graph, NULL, 2);
	vcn_graph_destroy(graph);

	double* F = calloc(2 * mesh->N_vertices, sizeof(*F));

	int status_assemble =
		pipeline_assemble_system(K, NULL, F, mesh, elemtype, material,
					 enable_self_weight, gravity,
					 analysis2D, params2D,
					 elements_enabled);
	if (0 != status_assemble)
		goto CLEANUP_LINEAR_SYSTEM;

	pipeline_set_boundary_conditions(mesh, K, F, bcond, 1.0);

  
	int solver_status = solver(K, F, displacement);
	if (0 != solver_status)
		goto CLEANUP_LINEAR_SYSTEM;

	pipeline_compute_strain(strain, mesh, displacement, elemtype,
				analysis2D, material);
	
	status = 0;
CLEANUP_LINEAR_SYSTEM:
	vcn_sparse_destroy(K);
	free(F);
	return status;
}

static inline int solver(const vcn_sparse_t *const A,
			 const double *const b, double* x)
{
	memset(x, 0, vcn_sparse_get_size(A) * sizeof(*x));
	return vcn_sparse_solve_CG_precond_Jacobi(A, b, x,
						  vcn_sparse_get_size(A),
						  1e-8, NULL, NULL, 1);
}

void vcn_fem_compute_stress_from_strain
			(uint32_t N_elements,
			 uint32_t* elements_connectivity_matrix, 
			 const vcn_fem_elem_t *const elemtype,
			 const vcn_fem_material_t *const material,
			 nb_analysis2D_t analysis2D,
			 double* strain,
			 const bool* elements_enabled /* NULL to enable all */,
			 double* stress /* Output */)
{
	/* Compute stress from element strain */
	uint32_t omp_parallel_threads = 1;
#pragma omp parallel for num_threads(omp_parallel_threads) schedule(guided)
	for(uint32_t i = 0; i < N_elements; i++){
		/* Get material properties */
		double E = vcn_fem_material_get_elasticity_module(material);

		if(elements_enabled != NULL)
			if(!elements_enabled[i])
				E = 0;

		double v = vcn_fem_material_get_poisson_module(material);
    
		/* Get constitutive matrix */
		double d11, d12, d22, d33;
		if (NB_PLANE_STRESS == analysis2D) {
			d11 = E / (1.0 - POW2(v));
			d12 = v * d11;
			d22 = d11;
			d33 = E / (2.0 * (1.0 + v));
		} else if (NB_PLANE_STRAIN == analysis2D) {
			d11 = (E * (1.0 - v)) / ((1.0 + v) * (1.0 - 2 * v));
			d12 = (v * d11) / (1.0 - v);
			d22 = d11;
			d33 = E / (2.0 * (1.0 + v));
		} else {
			/* Default: Plane stress */
			d11 = E / (1.0 - POW2(v));
			d12 = v * d11;
			d22 = d11;
			d33 = E / (2.0 * (1.0 + v));
		}
		/* Calculate stress */
		stress[i * 3] = strain[i * 3] * d11 + strain[i*3+1] * d12;
		stress[i*3+1] = strain[i * 3] * d12 + strain[i*3+1] * d22;
		stress[i*3+2] = strain[i*3+2] * d33;
	}
}

void vcn_fem_compute_von_mises(uint32_t N_elements,
			       double *stress,
			       double *von_mises /* Output */)
{
	for (uint32_t i = 0; i < N_elements; i++) {
		double sx = stress[i * 3];
		double sy = stress[i*3+1];
		double sxy = stress[i*3+2];
		von_mises[i] = sqrt(POW2(sx) + POW2(sy) - sx * sy +
				    3 * POW2(sxy));
	}
}
