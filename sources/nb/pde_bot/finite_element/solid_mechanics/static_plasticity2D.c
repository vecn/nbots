#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "nb/memory_bot.h"
#include "nb/math_bot.h"
#include "nb/io_bot.h"
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

//#include "nb/pde_bot/finite_element/solid_mechanics/set_bconditions.h"
#include "set_bconditions.h"
#include "pipeline.h"
//#include "nb/pde_bot/common_solid_mechanics/analysis2D.h"
#include "nb/pde_bot/finite_element/solid_mechanics/static_plasticity2D.h"
#include "nb/pde_bot/finite_element/solid_mechanics/plastic_stiffness_matrix.h"

#define POW2(a) ((a)*(a))

static void get_elem_regime(nb_plastified_analysis2D *elem_regime, uint32_t N_elem);

static void get_plastic_elements(uint32_t N_elem, nb_plastified_analysis2D *elem_regime, bool *plastic_elements);

/*void print_results_on_graph(uint32_t N_nod, uint32_t N_elem, double *total_displacement, double *stress,
                            nb_plastified_analysis2D *elem_regime, double *total_strain, const nb_mesh2D_t *const part,
                            bool *plastic_elements);*/

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
			 double *strain, /* Output */
			 double *stress, /*Output*/
			 double *displacement, /* Output */
			 uint32_t N_force_steps,
			 bool *plastic_elements, /* Output */
			 double *nodal_strain, /* Output */
			 double *nodal_stress /* Output */,
			 double *stiffness_factors,
			 double *density_factors)
{
	int status = 0;
	nb_graph_t *graph = nb_allocate_mem(nb_graph_get_memsize());
	nb_graph_init(graph);
	nb_mesh2D_load_graph(part, graph, NB_NODES_LINKED_BY_ELEMS);
	nb_sparse_t *K = nb_sparse_create(graph, NULL, 2);
	nb_graph_finish(graph);
	nb_free_mem(graph);

    	uint32_t N_elem = nb_mesh2D_get_N_elems(part);
	uint32_t N_nod = nb_mesh2D_get_N_nodes(part);

   	uint32_t F_elemsize = 3 * N_elem * sizeof(double);
    	uint32_t F_memsize = 2 * N_nod * sizeof(double);

    	uint32_t elem_regime_size = 10 * N_elem * sizeof(char);
    	uint32_t F_size = F_memsize;

    	uint64_t memsize = elem_regime_size + F_size;

    	char* memblock = malloc(memsize);

   	nb_plastified_analysis2D *elem_regime = (void*) memblock;
   	double *F = (void*)(memblock + elem_regime_size);
	memset(F, 0, F_memsize);
	memset(strain, 0, F_elemsize);
	memset(displacement, 0, F_memsize);

	get_elem_regime(elem_regime, N_elem);

    	double yield_stress = nb_material_get_yield_stress(material);

    	for (uint32_t i = 0; i < N_force_steps; i++) {
		double max_vm_stress = 0;
        	uint32_t N_plastic_elem = 0;
        	uint32_t plastified_elem = 0;

		uint8_t status_assemble = pipeline_assemble_plastic_system
								(K, NULL,F, part, elemtype, material,
		 						enable_self_weight, gravity, 
								analysis2D, params2D, elements_enabled, elem_regime, 
								stiffness_factors, density_factors);
		if (0 != status_assemble) {
			status = 1;
			goto CLEANUP_LINEAR_SYSTEM;
		}
	
		double condition_factor = (i + 1.0)/(double) N_force_steps;

		nb_fem_set_bconditions(part, K, F, bcond, condition_factor);
	
		
		double status = plastic_solver(K, F, displacement, N_nod);
	
		if(0 != status){
            		status = 2;
            		goto CLEANUP_LINEAR_SYSTEM;
        	}
        	printf("Iteration: %d\n", i); /* TEMPORAL */
		printf("Condition factor: %lf\n", condition_factor);
       		
        	pipeline_compute_strain(strain, part, displacement, elemtype);
		
        	nb_fem_compute_plastic_stress_from_strain(N_elem, elemtype, material, analysis2D,
                                                      strain, elements_enabled, stress, elem_regime);
		
        	get_stress_params(&max_vm_stress, &N_plastic_elem, &plastified_elem,
                          	stress, yield_stress, elem_regime, N_elem);
			
		while(N_plastic_elem != 0) {
			
        		max_vm_stress = nb_pde_get_vm_stress(stress[3*plastified_elem], stress[3*plastified_elem + 1], 
							stress[3*plastified_elem + 2]);
        		
			printf("Plastified elements: %d\n", N_plastic_elem);
			uint8_t status_assemble = pipeline_assemble_plastic_system
								(K, NULL,F, part, elemtype, material,
		 						enable_self_weight, gravity, 
								analysis2D, params2D, elements_enabled, elem_regime, 
								stiffness_factors, density_factors);
			if (0 != status_assemble) {
				status = 1;
				goto CLEANUP_LINEAR_SYSTEM;
			}

			nb_fem_set_bconditions(part, K, F, bcond, condition_factor);

			status = plastic_solver(K, F, displacement, N_nod);
			if(0 != status){
		    		status = 2;
		    		goto CLEANUP_LINEAR_SYSTEM;
			}

			pipeline_compute_strain(strain, part, displacement, elemtype);

			nb_fem_compute_plastic_stress_from_strain(N_elem, elemtype, material, analysis2D,
		                                            strain, elements_enabled, stress, elem_regime);

			find_new_plastic_elem(N_elem, elem_regime, stress, yield_stress, &max_vm_stress, &plastified_elem, &N_plastic_elem);
		}
    	}
	get_plastic_elements(N_elem, elem_regime, plastic_elements);
	nb_mesh2D_extrapolate_elems_to_nodes(part, 3, strain, nodal_strain);
	nb_mesh2D_extrapolate_elems_to_nodes(part, 3, stress, nodal_stress);

    	
    	//print_results_on_graph(N_nod, N_elem, displacement, stress, elem_regime, strain, part, plastic_elements);

   	CLEANUP_LINEAR_SYSTEM:
    	nb_free_mem(memblock);
    	nb_sparse_destroy(K);
		
	return status;
}

void nb_fem_compute_plastic_stress_from_strain
			(uint32_t N_elements,
			 const nb_fem_elem_t *const elem,
			 const nb_material_t *const material,
			 nb_analysis2D_t analysis2D,
			 double* strain,
			 const bool* elements_enabled,
			 double* stress ,
			 nb_plastified_analysis2D *elem_regime)
{

	/*uint32_t omp_parallel_threads = 1;
    	pragma omp parallel for num_threads(omp_parallel_threads) 		schedule(guided)*/
    	memset(stress, 0, 3*N_elements*sizeof(stress));
	for (uint32_t i = 0; i < N_elements; i++) {

		double De[4] = {1e-6, 1e-6, 1e-6, 1e-6};
		double Dp[4] = {1e-6, 1e-6, 1e-6, 1e-6};

        	if (pipeline_elem_is_enabled(elements_enabled, i)) {
        		switch(elem_regime[i]){
                		case NB_PLASTIC:{
                    			nb_pde_get_plastified_constitutive_matrix(Dp, material,
						       analysis2D, elem_regime[i]);
                    			nb_pde_get_constitutive_matrix(De, material, analysis2D);
                    			uint8_t N_gp = nb_fem_elem_get_N_gpoints(elem);
                    			for (int j = 0; j < N_gp; j++) {
                        			uint32_t id = i * N_gp + j;
                        			stress[id * 3] = strain[id * 3] * Dp[0] +
                            			strain[id*3 + 1] * Dp[1];
                        			stress[id*3 + 1] = strain[id * 3] * Dp[1] +
                            			strain[id*3 + 1] * Dp[2];
                       				stress[id*3 + 2] = strain[id*3 + 2] * Dp[3];
                    			}
                			break;
                		}
                		case NB_ELASTIC: {
                   			nb_pde_get_constitutive_matrix(De, material, analysis2D);
                   		 	uint8_t N_gp = nb_fem_elem_get_N_gpoints(elem);
                    			for (int j = 0; j < N_gp; j++) {
                      				uint32_t id = i * N_gp + j;
                        			stress[id * 3] = strain[id * 3] * De[0] +
                            			strain[id*3+1] * De[1];
                        			stress[id*3+1] = strain[id * 3] * De[1] +
                            			strain[id*3+1] * De[2];
                        			stress[id*3+2] = strain[id*3+2] * De[3];
                    			}
               				break;
                		}
            		}
        	}
	}
}

void get_stress_params(double *max_vm_stress, uint32_t *N_plastic_elem,
                       uint32_t *plastified_elem, double *stress, double yield_stress,
                       nb_plastified_analysis2D *elem_regime,
                       uint32_t N_elem) {
    	double vm_stress = 0;
    	max_vm_stress[0] = 0;
    	N_plastic_elem[0] = 1;
    	for(int j = 0; j < N_elem; j++) {
        	vm_stress = nb_pde_get_vm_stress(stress[3*j], stress[3*j +1], stress[3*j +2]);
        	if(vm_stress > yield_stress && elem_regime[j] != NB_PLASTIC) {
            		N_plastic_elem[0] += 1;
			elem_regime[j] = NB_PLASTIC;
            		if(vm_stress > max_vm_stress[0]) {
                		max_vm_stress[0] = vm_stress;
                		plastified_elem[0] = j;
            		}
        	}
    	}

    	if(max_vm_stress[0] == 0){
        	N_plastic_elem[0] = 0;
    	}
}


void find_new_plastic_elem(uint32_t N_elem,
                           nb_plastified_analysis2D *elem_regime,
                           double *stress, double yield_stress,
                           double *max_vm_stress,
                           uint32_t *plastified_elem,
                           uint32_t *N_plastic_elem)
{
    	plastified_elem[0] = 0;
	N_plastic_elem[0] = 0;

    	for(int j = 0; j < N_elem; j++) {
        	double vm_stress = nb_pde_get_vm_stress(stress[3*j], stress[3*j +1], stress[3*j +2]);
        	if(vm_stress > yield_stress && elem_regime[j] != NB_PLASTIC) {
            		N_plastic_elem[0] += 1;
			elem_regime[j] = NB_PLASTIC;
            		if(vm_stress > max_vm_stress[0] && elem_regime[j] != NB_PLASTIC) {
                		max_vm_stress[0] = vm_stress;
                		plastified_elem[0] = j;
            		}
        	}
    	}
}

int plastic_solver(const nb_sparse_t *const A,
		  const double *const b, double* x, uint32_t N_nod)
{
	uint32_t N = nb_sparse_get_size(A);
	memset(x, 0, 2 * N_nod * sizeof(double));
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

static void get_elem_regime(nb_plastified_analysis2D *elem_regime, uint32_t N_elem){
    	for (int i = 0; i < N_elem; i++) {
        	elem_regime[i] = NB_ELASTIC;
    	}
}

static void get_plastic_elements(uint32_t N_elem, nb_plastified_analysis2D *elem_regime, bool *plastic_elements) {
	for(int j = 0; j < N_elem; j++) {
        	switch(elem_regime[j]){
            		case NB_ELASTIC:
                		plastic_elements[j] = false;
                	break;
            		case NB_PLASTIC:
                		plastic_elements[j] = true;
               		break;
            		default:
                		plastic_elements[j] = false;
        	}
    	}
}

/*
void print_results_on_graph(uint32_t N_nod, uint32_t N_elem, double *total_displacement, double *stress,
                       nb_plastified_analysis2D *elem_regime, double *total_strain, const nb_mesh2D_t *const part,
                       bool *plastic_elements)
{
    uint32_t nodesize = N_nod * sizeof(double);
    uint32_t elemsize = N_elem * sizeof(double);

    uint32_t avg_displacement_size = nodesize;
    uint32_t total_displacement_x_size = nodesize;
    uint32_t total_displacement_y_size = nodesize;
    uint32_t Sx_size = elemsize;
    uint32_t Sy_size = elemsize;
    uint32_t Sxy_size = elemsize;
    //uint32_t plastic_elements_size = N_elem * sizeof(bool);

    uint64_t memsize = avg_displacement_size + total_displacement_x_size + total_displacement_y_size +
                        Sx_size + Sy_size + Sxy_size;// + plastic_elements_size;

    char* memblock = malloc(memsize);

    double *avg_displacement = (void*)memblock; //nb_allocate_mem(nodesize);
    for(int i = 1; i < N_nod; i++){
        avg_displacement[i] = sqrt(POW2(total_displacement[2*i])+ POW2(total_displacement[2*i +1]));
    }

    double *total_displacement_x = (void*)(memblock + avg_displacement_size); //nb_allocate_mem(nodesize);
    double *total_displacement_y = (void*)(memblock + avg_displacement_size + total_displacement_x_size); //nb_allocate_mem(nodesize);

    for(int j = 0; j < N_nod; j++){
        total_displacement_x[j] = total_displacement[2*j];
        total_displacement_y[j] = total_displacement[2*j +1];
    }
    double *Sx = (void*)(memblock + avg_displacement_size + total_displacement_x_size + total_displacement_y_size); //nb_allocate_mem(elemsize);
    double *Sy =(void*)(memblock + avg_displacement_size + total_displacement_x_size + total_displacement_y_size + Sx_size); // nb_allocate_mem(elemsize);
    double *Sxy =(void*)(memblock + avg_displacement_size + total_displacement_x_size + total_displacement_y_size + Sx_size +
                          Sy_size); // nb_allocate_mem(elemsize);
    for(int j = 0; j < N_elem; j++) {
        Sx[j] = stress[3*j];
        Sy[j] = stress[3*j + 1];
        Sxy[j] = stress[3*j +2];
    }
    
    for(int j = 0; j < N_elem; j++) {
        switch(elem_regime[j]){
            case NB_ELASTIC:
                plastic_elements[j] = false;
                break;
            case NB_PLASTIC:
                plastic_elements[j] = true;
                break;
            default:
                plastic_elements[j] = false;
        }
    }
    //Position of the palette is controlled by float label_width in static void add_palette() of draw.c
    //nb_mesh2D_distort_with_field(part, NB_NODE, total_displacement, 0.5);
    nb_mesh2D_export_draw(part, "cantilever_beam_plastic_elements.png", 1200, 800, NB_ELEMENT, NB_CLASS, plastic_elements, true);
    nb_mesh2D_export_draw(part,"cantilever_beam_avg_disp.png", 1200, 800, NB_NODE, NB_FIELD, avg_displacement, true);
    nb_mesh2D_export_draw(part,"cantilever_beam_disp_x.png", 1200, 800, NB_NODE, NB_FIELD, total_displacement_x, true);
    nb_mesh2D_export_draw(part,"cantilever_disp_y.png", 1200, 800, NB_NODE, NB_FIELD, total_displacement_y, true);
    nb_mesh2D_export_draw(part,"cantilever_Sx.png", 1200, 800, NB_ELEMENT, NB_FIELD, Sx, true);
    nb_mesh2D_export_draw(part,"cantilever_Sy.png", 1200, 800, NB_ELEMENT, NB_FIELD, Sy, true);
    nb_mesh2D_export_draw(part,"cantilever_Sxy.png", 1200, 800, NB_ELEMENT, NB_FIELD, Sxy, true);

    nb_free_mem(memblock);
    //nb_free_mem(avg_displacement);
    //nb_free_mem(total_displacement_x);
    //nb_free_mem(total_displacement_y);
    //nb_free_mem(Sx);
    //nb_free_mem(Sy);
    //nb_free_mem(Sxy);
    //nb_free_mem(plastic_elements);
}*/
void interpolate_trg_strain_nodes_to_elemes(const nb_mesh2D_t *const part, double *strain, double *nodal_strain) 
{
	uint32_t N_elems = nb_mesh2D_get_N_elems(part);
	memset(strain, 0, 3 * N_elems * sizeof(*strain));

	for(int i = 0; i < N_elems; i++) {
		strain[3 * i] = nodal_strain[3 * nb_mesh2D_elem_get_adj(part, i, 0)]/3 + 
					nodal_strain[3 * nb_mesh2D_elem_get_adj(part, i, 1)]/3 + 
					nodal_strain[3 * nb_mesh2D_elem_get_adj(part, i, 2)]/3;

		strain[3 * i + 1] = nodal_strain[3 * nb_mesh2D_elem_get_adj(part, i, 0) + 1]/3 + 
					nodal_strain[3 * nb_mesh2D_elem_get_adj(part, i, 1) + 1]/3 + 
					nodal_strain[3 * nb_mesh2D_elem_get_adj(part, i, 2) + 1]/3;

		strain[3 * i + 2] = nodal_strain[3 * nb_mesh2D_elem_get_adj(part, i, 0) + 2]/3 + 
					nodal_strain[3 * nb_mesh2D_elem_get_adj(part, i, 1) + 2]/3 + 
					nodal_strain[3 * nb_mesh2D_elem_get_adj(part, i, 2) + 2]/3;	
	}
}

void interpolate_trg_stress_nodes_to_elemes(const nb_mesh2D_t *const part, double *stress, double *nodal_stress)
{
	uint32_t N_elems = nb_mesh2D_get_N_elems(part);
	memset(stress, 0, 3 * N_elems * sizeof(*stress));

	for(int i = 0; i < N_elems; i++) {
		stress[3 * i] = nodal_stress[3 * nb_mesh2D_elem_get_adj(part, i, 0)]/3 + 
					nodal_stress[3 * nb_mesh2D_elem_get_adj(part, i, 1)]/3 + 
					nodal_stress[3 * nb_mesh2D_elem_get_adj(part, i, 2)]/3;

		stress[3 * i + 1] = nodal_stress[3 * nb_mesh2D_elem_get_adj(part, i, 0) + 1]/3 + 
					nodal_stress[3 * nb_mesh2D_elem_get_adj(part, i, 1) + 1]/3 + 
					nodal_stress[3 * nb_mesh2D_elem_get_adj(part, i, 2) + 1]/3;

		stress[3 * i + 2] = nodal_stress[3 * nb_mesh2D_elem_get_adj(part, i, 0) + 2]/3 + 
					nodal_stress[3 * nb_mesh2D_elem_get_adj(part, i, 1) + 2]/3 + 
					nodal_stress[3 * nb_mesh2D_elem_get_adj(part, i, 2) + 2]/3;	
	}
}
