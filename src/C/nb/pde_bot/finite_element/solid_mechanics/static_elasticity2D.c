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
#include "nb/pde_bot/boundary_conditions.h"
#include "nb/pde_bot/finite_element/element.h"
#include "nb/pde_bot/finite_element/gaussp_to_nodes.h"
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
			 const vcn_bcond_t *const bcond,
			 char enable_self_weight,
			 double gravity[2],
			 bool enable_plane_stress,
			 double thickness,
			 const bool *elements_enabled, /* NULL to enable all */
			 double *displacement, /* Output */
			 double *strain,       /* Output */
			 const char* logfile    /* NULL if not required */)
/* The output vectors must be allocated before start:
 *     > displacement:  2 * N_vertices (size of double)
 *     >       strain:  3 * N_vertices (size of double)
 */
{
	int status = 1;
	if (NULL != logfile) {
		FILE *log = fopen(logfile, "w");
		fprintf(log, "FINITE ELEMENT\n");
		fprintf(log, "    > SOLID MECHANICS\n");
		fprintf(log, "        > 2D STATIC ELASTICITY\n\n");
		fclose(log);
	}

	/*********************************************************************/
	/****************** 1) Assemble system *******************************/
	/*********************************************************************/
	/* Allocate global Stiffness Matrix and Force Vector */
	vcn_graph_t *graph = vcn_msh3trg_create_vtx_graph(mesh);
	vcn_sparse_t *K = vcn_sparse_create(graph, NULL, 2);
	vcn_graph_destroy(graph);

	double* F = calloc(2 * mesh->N_vertices, sizeof(*F));
	/* Allocate elemental Stiffness Matrix and Force Vector */
	int status_assemble =
		pipeline_assemble_system(K, NULL, F, mesh, elemtype, material,
					 enable_self_weight, gravity,
					 enable_plane_stress, thickness,
					 elements_enabled);
	if (0 != status_assemble) {
		if (NULL != logfile) {
			FILE* log = fopen(logfile, "a");
			fprintf(log, "Assemble system fails (Code: %i).\n",
				status_assemble);
			fclose(log);
		}
		goto CLEANUP_LINEAR_SYSTEM;
	}

	/*********************************************************************/
	/**************** 2) Set boundary conditions *************************/
	/*********************************************************************/
	pipeline_set_boundary_conditions(mesh, K, F, bcond, thickness, 1.0);

	/**********************************************************************/
	/************* 3) Solve system (to compute displacements) *************/
	/**********************************************************************/
  
	int solver_status = solver(K, F, displacement);
	/* Display failure info in logfile */
	if (0 != solver_status) {
		if (NULL != logfile) {
			FILE* log = fopen(logfile, "a");
			fprintf(log, "Solver fails (Code: %i).\n", solver_status);
			fclose(log);
		}
		goto CLEANUP_LINEAR_SYSTEM;
	}

	/********************************************************************/
	/********************* 4) Compute Strain            *****************/
	/*********************       >  S = B u             *****************/
	/********************************************************************/
	pipeline_compute_strain(strain, mesh, displacement, elemtype,
				enable_plane_stress, material);
	
	status = 0;
CLEANUP_LINEAR_SYSTEM:
	vcn_sparse_destroy(K);
	free(F);
	return status;
}

static inline int solver(const vcn_sparse_t *const A,
			 const double *const b, double* x)
{
	return vcn_sparse_solve_CG_precond_Jacobi(A, b, x,
						  vcn_sparse_get_size(A),
						  1e-8, NULL, NULL, 1);
}

void vcn_fem_compute_stress_from_strain
			(uint32_t N_elements,
			 uint32_t* elements_connectivity_matrix, 
			 const vcn_fem_elem_t *const elemtype,
			 const vcn_fem_material_t *const material,
			 bool enable_plane_stress,
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
		double d11 = E/(1.0 - POW2(v));
		double d12 = v*d11;
    
		if(!enable_plane_stress){
			d11 = (E*(1.0-v))/((1.0 + v)*(1.0-2*v));
			d12 = (v*d11)/(1.0-v);
		}    
		double d22 = d11;
		double d33 = E/(2.0*(1.0+v));
		/* Calculate stress */
		stress[i * 3] = strain[i * 3] * d11 + strain[i*3+1] * d12;
		stress[i*3+1] = strain[i * 3] * d12 + strain[i*3+1] * d22;
		stress[i*3+2] = strain[i*3+2] * d33;
	}
}
