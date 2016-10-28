#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "nb/memory_bot.h"
#include "nb/math_bot.h"
#include "nb/cfreader_cat.h"
#include "nb/solver_bot.h"
#include "nb/container_bot.h"
#include "nb/graph_bot.h"
#include "nb/geometric_bot.h"
#include "nb/pde_bot/material.h"
#include "nb/pde_bot/common_solid_mechanics/analysis2D.h"
#include "nb/pde_bot/common_solid_mechanics/formulas.h"
#include "nb/pde_bot/boundary_conditions/bcond.h"
#include "nb/pde_bot/boundary_conditions/bcond_iter.h"
#include "nb/pde_bot/finite_element/element.h"
#include "nb/pde_bot/finite_element/gaussp_to_nodes.h"
#include "nb/pde_bot/finite_element/solid_mechanics/static_elasticity2D.h"

#include "nb/pde_bot/finite_element/solid_mechanics/set_bconditions.h"
#include "pipeline.h"
#include "nb/pde_bot/common_solid_mechanics/analysis2D.h"
#include "nb/pde_bot/finite_element/solid_mechanics/static_plasticity2D.h"

#define POW2(a) ((a)*(a))

double stress_max_tolerance(double *plastic_stress, nb_plastified_analysis2D *elem_regime, uint32_t N_elem);

int fem_compute_plastic_2D_Solid_Mechanics
			(const nb_mesh2D_t *const part,
			 const nb_fem_elem_t *const elemtype,
			 const nb_material_t *const material,
			 const nb_bcond_t *const bcond,
			 bool enable_self_weight,
			 double gravity[2],
			 nb_analysis2D_t analysis2D,
			 nb_analysis2D_params *params2D,
			 const bool *elements_enabled, /* NULL to enable all */
			 double *total_strain, /* Output*/
			 double *stress,
			 double *displacement, /* Output, just the last computed plastic displacement */
			 uint32_t N_force_steps,
			 double accepted_tol)
{
    printf("Number of force steps: %i\n", N_force_steps);
    printf("Accepted tolerance for the equilibrium of forces: %f\n", accepted_tol);
	int status = 0;
	nb_graph_t *graph = malloc(nb_graph_get_memsize());
	nb_graph_init(graph);
	nb_mesh2D_load_graph(part, graph, NB_NODES_LINKED_BY_ELEMS);
	nb_sparse_t *K = nb_sparse_create(graph, NULL, 2);
	nb_graph_finish(graph);

	uint32_t N_nod = nb_mesh2D_get_N_nodes(part);
	printf("Number of nodes: %i\n", N_nod);
	uint32_t N_elem = nb_mesh2D_get_N_elems(part);
	printf("Number of elements: %i\n", N_elem);
	uint32_t F_memsize = 2 * N_nod * sizeof(double);
	uint32_t F_elemsize = 3*N_elem*sizeof(double);
	double *F = nb_soft_allocate_mem(F_memsize);
	memset(F, 0, F_memsize);
	memset(total_strain, 0, F_elemsize);
	double *el_strain = nb_allocate_zero_mem(F_elemsize);
	double *strain = nb_allocate_zero_mem(F_elemsize);
	double *elastic_strain = nb_allocate_zero_mem(F_elemsize);
	double *last_elastic_strain = nb_allocate_zero_mem(F_elemsize);
    double *plastic_strain = nb_allocate_zero_mem(F_elemsize);
    double *elastic_stress = nb_allocate_zero_mem(F_elemsize);
    double *plastic_stress = nb_allocate_zero_mem(F_elemsize);

	int status_assemble =
	pipeline_assemble_system(K, NULL, F, part, elemtype, material,
					 enable_self_weight, gravity,
					 analysis2D, params2D,
					 elements_enabled);
	if (0 != status_assemble) {
		status = 1;
		goto CLEANUP_LINEAR_SYSTEM;
	}

	nb_fem_set_bconditions(part, K, F, bcond, 1.0);

    nb_plastified_analysis2D *elem_regime = nb_allocate_mem(N_elem*10*sizeof(char));

    for (int i = 0; i < N_elem; i++) {
        elem_regime[i] = NB_ELASTIC;
    }

    bool regime_status = false;//false = Elastic regime, true = Plastic regime

    double *dF_basic = nb_soft_allocate_mem(F_memsize);
    memset(dF_basic, 0, F_memsize);
	for (int i = 0; i < 2*N_nod; i++) {
        dF_basic[i] = F[i] / N_force_steps;
    }
    double *dF = nb_soft_allocate_mem(F_memsize);
    memset(dF, 0, F_memsize);


    /* Apply the force recursively step by step */
    for (int i = 0; i < N_force_steps; i++) {
        for (int j = 0; j < 2*N_nod; j++)
        {
            dF[j] += dF_basic[j];
        }
    /* TEMPORAL
    for(int k = 0; k < N_nod; k++) {
        printf("dFx[%d]: %lf \t", k, dF[2*k]);
        printf("dFy[%d]: %lf \n", k, dF[2*k+1]);
    }*/

    int solver_status = plastic_solver(K, dF, displacement);
    /* TEMPORAL
    for(int k= 0; k < N_nod; k++) {
        printf("Dx[%d]: %lf \t ", k, displacement[2*k]);
        printf("Dy[%d]: %lf \n ", k, displacement[2*k + 1]);
    }*/

    if (0 != solver_status) {
        status = 2;
        goto CLEANUP_LINEAR_SYSTEM;
    }

    /* Compute strains of the elastic displacements. The strains are stored in the elastic strain vector if and just if the
    element has not plastified yet: */
	add_elastic_strain(el_strain, last_elastic_strain, elem_regime, N_elem);

	pipeline_compute_strain(el_strain, part, displacement, elemtype);
    add_elastic_strain(el_strain, elastic_strain, elem_regime, N_elem);
    add_total_strain(total_strain, elastic_strain, plastic_strain, N_elem);

	/* Compute stresses from total strain. The elastic elements stresses are computed with the elastic constitutive matrix and
	the plastic elements with the plastic constitutive matrix with the plastic module Ep: */
	nb_fem_compute_plastic_stress_from_strain(N_elem, elemtype, material, analysis2D, total_strain,
                                              elements_enabled, stress, elem_regime);

    /* Compute Von Mises stress for each non plastified element:*/
    double yield_stress = nb_material_get_yield_stress(material);
    double vm_stress;
    for (int j = 0; j < N_elem; j++) {
        if (elem_regime[j] != NB_PLASTIC) {
            vm_stress = nb_pde_get_vm_stress(stress[3*j], stress[3*j+1], stress[3*j+2]);
            if (vm_stress > yield_stress) {
                elem_regime[j] = NB_PLASTIC;
                if(!regime_status) {
                    regime_status = true;
                }
            }
        }
    }

    double *FI = nb_soft_allocate_mem(F_memsize);

    if (regime_status) {
        nb_compute_plastic_strain(strain, total_strain, last_elastic_strain, N_elem, elem_regime);
        /* Recalculate the stresses of the elements taking into account the new plastified elements: */
        nb_fem_compute_plastic_stress_from_strain(N_elem, elemtype, material, analysis2D, strain,
                                                elements_enabled, plastic_stress, elem_regime);
        nb_fem_compute_elastic_stress_from_strain(N_elem, elemtype, material, analysis2D, strain,
                                                elements_enabled ,elastic_stress);

        double *res_stress = nb_allocate_mem(F_elemsize); /* LIBERAR MEMORIA*/

        nb_fem_compute_diference_of_stresses(res_stress, elastic_stress, plastic_stress, N_elem, elem_regime);
        /* TEMPORAL*/
        for (int j = 0; j < N_elem; j++) {
                vm_stress = nb_pde_get_vm_stress(res_stress[3*j], res_stress[3*j+1], res_stress[3*j+2]);
                printf("Von Misses [%d]: %lf\n", j, vm_stress);
        }
        /* TEMPORAL */
        for(int k = 0; k < N_elem; k++) {
            printf("Sx[%d]: %lf\t\t", k, res_stress[3*k]);
            printf("Sy[%d]: %lf\t\t", k, res_stress[3*k+1]);
            printf("Sxy[%d]: %lf\n", k, res_stress[3*k+2]);
        }

        /* Compute internal forces: */
        memset(FI, 0, F_memsize);
        int status_internal = compute_internal_forces(FI, res_stress, N_elem, N_nod, part,
                                                      material, analysis2D, params2D, elemtype);
        /*TEMPORAL*/
        for (int k = 0; k < N_nod; k++) {
            printf("FIx[%d]: %lf\t", k, FI[2*k]);
            printf("FIy[%d]: %lf\n", k, FI[2*k+1]);
        }
        if(0 != status_internal) {
            status = 3;
            goto CLEANUP_RECURSIVE_SYSTEM;
        }

        double tolerance = stress_max_tolerance(res_stress, elem_regime, N_elem);

        while(tolerance > accepted_tol) {
            int solver_status = plastic_solver(K, FI, displacement);
            if (0 != solver_status) {
                status = 2;
                goto CLEANUP_LINEAR_SYSTEM;
            }
            pipeline_compute_strain(strain, part, displacement, elemtype);
            add_plastic_strain(strain, plastic_strain, elem_regime, N_elem);
            add_total_strain(total_strain, elastic_strain, plastic_strain, N_elem);

            nb_fem_compute_plastic_stress_from_strain(N_elem, elemtype, material, analysis2D, strain,
                                                elements_enabled, plastic_stress, elem_regime);

            nb_fem_compute_diference_of_stresses(res_stress, elastic_stress, plastic_stress, N_elem, elem_regime);

                    /* TEMPORAL*/
            for (int j = 0; j < N_elem; j++) {
                    vm_stress = nb_pde_get_vm_stress(res_stress[3*j], res_stress[3*j+1], res_stress[3*j+2]);
                    printf("Von Misses [%d]: %lf\n", j, vm_stress);
                    break;
            }
        memset(FI, 0, F_memsize);
        int status_internal = compute_internal_forces(FI, res_stress, N_elem, N_nod, part,
                                                      material, analysis2D, params2D, elemtype);
        if(0 != status_internal) {
            status = 3;
            goto CLEANUP_RECURSIVE_SYSTEM;
        }

        double tolerance = stress_max_tolerance(res_stress, elem_regime, N_elem);
        }
        /* Compute the maximum difference between the internal and the external forces as the tolerance: */
       // double tolerance = force_tolerance(dF, FI, N_nod);
        /* Recursive process to balance the actual external forces with internal ones */
        while(tolerance > accepted_tol) {
            double *res_force = nb_soft_allocate_mem(F_memsize);
            memset(res_force, 0, F_memsize);
            /* Compute the residual force vector between external dF forces and internal FI forces*/
            residual_force(res_force, dF, FI, N_nod);
            /* Solve the linear system of equations with the residual forces and accumulate
            the displacements */
            int solver_status = plastic_solver(K, FI, displacement);

            nb_soft_free_mem(F_memsize, res_force);
            if (0 != solver_status) {
                status = 2;
                goto CLEANUP_LINEAR_SYSTEM;
            }
            /* Compute the strains from the computed displacements and add them to the plastic strains vectors and to the
            global strain vector afterwards: */
            pipeline_compute_strain(strain, part, displacement, elemtype);
            add_plastic_strain(strain, plastic_strain, elem_regime, N_elem);
            add_total_strain(total_strain, elastic_strain, plastic_strain, N_elem);

            nb_fem_compute_plastic_stress_from_strain(N_elem, elemtype, material, analysis2D, total_strain,
                                                            elements_enabled, stress, elem_regime);

            /* Recompute internal forces from stresses and setup the new tolerance: */
            memset(FI, 0, F_memsize);
            int status_internal = compute_internal_forces(FI, stress, N_elem, N_nod, part,
                                                          material, analysis2D, params2D, elemtype);
            if(0 != status_internal) {
                status = 3;
                goto CLEANUP_RECURSIVE_SYSTEM;
            }
            //double tolerance = force_tolerance(dF, FI, N_nod);
        }

    }
    CLEANUP_RECURSIVE_SYSTEM:
    nb_soft_free_mem(F_memsize, FI);
    if (status == 3){
    goto CLEANUP_LINEAR_SYSTEM;
    }
}

CLEANUP_LINEAR_SYSTEM:
	nb_sparse_destroy(K);
	nb_soft_free_mem(F_memsize, F);
	nb_soft_free_mem(F_memsize, dF_basic);
	nb_soft_free_mem(F_memsize, dF);
	nb_free_mem(elastic_strain);
	nb_free_mem(plastic_strain);
	nb_free_mem(strain);
	nb_free_mem(elem_regime);

	return status;
}
void nb_fem_compute_plastic_stress_from_strain
			(uint32_t N_elements,
			 const nb_fem_elem_t *const elem,
			 const nb_material_t *const material,
			 nb_analysis2D_t analysis2D,
			 double* strain,
			 const bool* elements_enabled ,
			 double* stress ,
			 nb_plastified_analysis2D *elem_regime)
{

	//uint32_t omp_parallel_threads = 1;
    #pragma omp parallel for num_threads(omp_parallel_threads) schedule(guided)
    memset(stress, 0, 3*N_elements*sizeof(stress));
	for (uint32_t i = 0; i < N_elements; i++) {

		double D[4] = {1e-6, 1e-6, 1e-6, 1e-6};

		if (pipeline_elem_is_enabled(elements_enabled, i))
        nb_pde_get_plastified_constitutive_matrix(D, material,
						       analysis2D, elem_regime[i]);
		uint8_t N_gp = nb_fem_elem_get_N_gpoints(elem);
		for (int j = 0; j < N_gp; j++) {
			uint32_t id = i * N_gp + j;
			stress[id * 3] = strain[id * 3] * D[0] +
				strain[id*3+1] * D[1];
			stress[id*3+1] = strain[id * 3] * D[1] +
				strain[id*3+1] * D[2];
			stress[id*3+2] = strain[id*3+2] * D[3];
		}
	}
}

void nb_fem_compute_elastic_stress_from_strain
			(uint32_t N_elements,
			 const nb_fem_elem_t *const elem,
			 const nb_material_t *const material,
			 nb_analysis2D_t analysis2D,
			 double* strain,
			 const bool* elements_enabled  ,
			 double* stress )
{
	//uint32_t omp_parallel_threads = 1;
    #pragma omp parallel for num_threads(omp_parallel_threads) schedule(guided)
    memset(stress, 0, 3*N_elements*sizeof(stress));
	for (uint32_t i = 0; i < N_elements; i++) {

		double D[4] = {1e-6, 1e-6, 1e-6, 1e-6};

		if (pipeline_elem_is_enabled(elements_enabled, i)) {
            nb_pde_get_constitutive_matrix(D, material, analysis2D);
            uint8_t N_gp = nb_fem_elem_get_N_gpoints(elem);
            for (int j = 0; j < N_gp; j++) {
                uint32_t id = i * N_gp + j;
                stress[id * 3] = strain[id * 3] * D[0] +
                    strain[id*3+1] * D[1];
                stress[id*3+1] = strain[id * 3] * D[1] +
                    strain[id*3+1] * D[2];
                stress[id*3+2] = strain[id*3+2] * D[3];
            }
		}
    }
}

void nb_fem_compute_diference_of_stresses(double *residual_stress, double *elastic_stress,
                                          double *plastic_stress, uint32_t N_elements,
                                          nb_plastified_analysis2D *elem_regime) {
    memset(residual_stress, 0, 3*N_elements*sizeof(residual_stress));
    for(int i = 0; i < N_elements; i++) {
        switch(elem_regime[i]){
        case NB_ELASTIC:
            residual_stress[3*i] = 0;
            residual_stress[3*i + 1] = 0;
            residual_stress[3*i + 2] = 0;
        case NB_PLASTIC:
            residual_stress[3*i] = elastic_stress[3*i] - plastic_stress[3*i];
            residual_stress[3*i + 1] = elastic_stress[3*i + 1] - plastic_stress[3*i + 1];
            residual_stress[3*i + 2] = elastic_stress[3*i + 2] - plastic_stress[3*i + 2];
        }
    }
}

void residual_force(double *res_force, double *dF, double *FI, uint32_t N_nod) {

    for (int i = 0; i < 2*N_nod; i++) {
        res_force[i] = dF[i] - FI[i];
    }
}

void add_elastic_strain(double *strain, double *elastic_strain, nb_plastified_analysis2D *elem_regime, uint32_t N_elem) {

    for (int i = 0; i < N_elem; i++) {
        switch(elem_regime[i]){
        case NB_ELASTIC:
            elastic_strain[3*i] = strain[3*i];
            elastic_strain[3*i+1] = strain[3*i+1];
            elastic_strain[3*i+2] = strain[3*i+2];
            break;
        case NB_PLASTIC:
            elastic_strain[3*i] += 0;
            elastic_strain[3*i+1] += 0;
            elastic_strain[3*i+2] += 0;
            break;
        }
    }
}

void add_plastic_strain(double *strain, double *plastic_strain, nb_plastified_analysis2D *elem_regime, uint32_t N_elem) {

    for (int i = 0; i < N_elem; i += 3) {
        switch(elem_regime[i]){
        case NB_ELASTIC:
            plastic_strain[3*i] += 0;
            plastic_strain[3*i + 1] += 0;
            plastic_strain[3*i + 2] += 0;
            break;
        case NB_PLASTIC:
            plastic_strain[3*i] += strain[3*i];
            plastic_strain[3*i + 1] += strain[3*i + 1];
            plastic_strain[3*i + 2] += strain[3*i + 2];
            break;
        }
    }
}

void add_total_strain(double *total_strain, double *elastic_strain, double *plastic_strain, uint32_t N_elem) {

    for (int i = 0; i < 3*N_elem; i++) {
        total_strain[i] = elastic_strain[i] + plastic_strain[i];
    }
}

int plastic_solver(const nb_sparse_t *const A,
		  const double *const b, double* x)
{
	uint32_t N = nb_sparse_get_size(A);
	memset(x, 0, N * sizeof(*x));
	int status = nb_sparse_solve_CG_precond_Jacobi(A, b, x, N,
							1e-8, NULL,
							NULL, 1);
	int out;
	if (0 == status || 1 == status)
		out = 0;
	else
		out = 1; // Tolerance not reached in CG Jacobi
	return out;
}

void nb_compute_plastic_strain(double *plastic_strain, double *total_strain, double *last_elastic_strain,
                               uint32_t N_elem, nb_plastified_analysis2D *elem_regime) {
    for (int i = 0; i < N_elem; i++) {
        switch(elem_regime[i]){
        case NB_ELASTIC:
        plastic_strain[3*i] = 0;
        plastic_strain[3*i + 1] = 0;
        plastic_strain[3*i + 2] = 0;
        case NB_PLASTIC:
        plastic_strain[3*i] = total_strain[3*i] - last_elastic_strain [3*i];
        plastic_strain[3*i + 1] = total_strain[3*i + 1] - last_elastic_strain [3*i + 1];
        plastic_strain[3*i + 2] = total_strain[3*i + 2] - last_elastic_strain [3*i + 2];
        break;
        }
    }
}

double stress_max_tolerance(double *plastic_stress, nb_plastified_analysis2D *elem_regime, uint32_t N_elem) {
    double vm_stress = 0;
    double tolerance = 0;
    for (int j = 0; j < N_elem; j++) {
        switch(elem_regime[j]){
        case NB_ELASTIC:
            vm_stress = 0;
            break;
        case NB_PLASTIC:
            vm_stress = nb_pde_get_vm_stress(plastic_stress[3*j], plastic_stress[3*j + 1], plastic_stress[3*j +2]);
            if(vm_stress > tolerance){
                tolerance = vm_stress;
            }
            break;
        }
    }
    return tolerance;
}

