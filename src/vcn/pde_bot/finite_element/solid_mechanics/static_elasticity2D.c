#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "vcn/math_bot.h"
#include "vcn/cfreader_cat.h"
#include "vcn/eigen_bot.h"
#include "vcn/container_bot.h"
#include "vcn/graph_bot.h"
#include "vcn/pde_bot/material.h"
#include "vcn/pde_bot/boundary_conditions.h"
#include "vcn/pde_bot/finite_element/element.h"
#include "vcn/pde_bot/finite_element/gaussp_to_nodes.h"
#include "vcn/pde_bot/finite_element/solid_mechanics/static_elasticity2D.h"
  
#include "../element_struct.h"
#include "pipeline.h"

#define POW2(a) ((a)*(a))

int vcn_fem_compute_2D_Solid_Mechanics
			(const vcn_msh3trg_t *const mesh,
			 const vcn_fem_elem_t *const elemtype,
			 const vcn_fem_material_t *const material,
			 const vcn_bcond_t *const bmeshcond,
			 char enable_self_weight,
			 double gravity[2],
			 int (*solver)(const vcn_sparse_t *const A,
				       const double *const b,
				       double* x, uint32_t omp_threads),
			 bool enable_plane_stress,
			 double thickness,
			 uint32_t omp_parallel_threads,
			 bool* elements_enabled, /* NULL to enable all */
			 double * displacement, /* Output */
			 double * strain,       /* Output */
			 const char* logfile    /* NULL if not required */)
/* The output vectors must be allocated before start:
 *     > displacement:  2 * N_vertices (size of double)
 *     >       strain:  3 * N_vertices (size of double)
 */
{
	if (NULL != logfile) {
		FILE *log = fopen(logfile, "a");
		fprintf(log, "Bidimensional elastic Solid Mechanics Solver using FEM\n");
		fclose(log);
	}

	/****************************************************************************/
	/********************** 1) Assemble system **********************************/
	/****************************************************************************/
	/* Allocate global Stiffness Matrix and Force Vector */
	vcn_graph_t *graph = vcn_msh3trg_create_vtx_graph(mesh);
	vcn_sparse_t *K = vcn_sparse_create(graph, NULL, 2);
	vcn_graph_destroy(graph);

	double* F = (double*)calloc(2 * mesh->N_vertices, sizeof(double));
	/* Allocate elemental Stiffness Matrix and Force Vector */
	pipeline_assemble_system(K, NULL, F,
				 mesh,
				 elemtype,
				 material,
				 enable_self_weight,
				 gravity,
				 enable_plane_stress,
				 thickness,
				 elements_enabled);
	/*****************************************************************************/
	/********************** 2) Set boundary conditions ***************************/
	/*****************************************************************************/
	pipeline_set_boundary_conditions(K, F, bmeshcond, thickness, 1.0);

	/*****************************************************************************/
	/**************** 3) Solve system (to compute displacements) *****************/
	/*****************************************************************************/
  
	char solver_status = solver(K, F, displacement, omp_parallel_threads);
    
	/* Display failure info in logfile */
	if(solver_status != 0){
		if(logfile != NULL){
			FILE* log = fopen(logfile, "a");
			fprintf(log, "Solver fails (Code: %i).\n", solver_status);
			fclose(log);
		}
		vcn_sparse_destroy(K);
		free(F);
		return 1;
	}

	vcn_sparse_destroy(K);
	free(F);
	/*****************************************************************************/
	/********************* 4) Compute Strain            **************************/
	/*********************       >  S = B u             **************************/
	/*****************************************************************************/
	pipeline_compute_strain(strain, mesh, displacement, elemtype,
				enable_plane_stress, material);
	/* Successful exit */
	return 0;
}

void vcn_fem_compute_stress_from_strain
(uint32_t N_elements,
 uint32_t* elements_connectivity_matrix, 
 const vcn_fem_elem_t *const elemtype,
 const vcn_fem_material_t *const material,
 bool enable_plane_stress,
 double* strain,
 uint32_t omp_parallel_threads,
 bool* elements_enabled /* NULL to enable all */,
 double* stress /* Output */){
	/* Compute stress from element strain */
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
