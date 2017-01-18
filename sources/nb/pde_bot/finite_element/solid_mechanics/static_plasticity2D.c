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

//#include "nb/pde_bot/finite_element/solid_mechanics/set_bconditions.h"
#include "set_bconditions.h"
#include "pipeline.h"
//#include "nb/pde_bot/common_solid_mechanics/analysis2D.h"
#include "nb/pde_bot/finite_element/solid_mechanics/static_plasticity2D.h"
#include "nb/pde_bot/finite_element/solid_mechanics/plastic_stiffness_matrix.h"

#define POW2(a) ((a)*(a))

static void get_elem_regime(nb_plastified_analysis2D *elem_regime, uint32_t N_elem);

static void get_dF_basic(double *dF_basic, uint32_t F_memsize, uint32_t N_nod, double *F, uint32_t N_force_steps);

static int get_first_plastic_element(const nb_sparse_t *const K, double *dF_basic, double *displacement, uint32_t N_nod,
                              double *total_displacement, const nb_mesh2D_t *const part, double *strain,
                              const nb_fem_elem_t *const elemtype, const nb_material_t *const material,
                              nb_analysis2D_t analysis2D, const bool *elements_enabled, double *stress,
                              nb_plastified_analysis2D *elem_regime, double *elastic_strain,
                              uint32_t i, uint32_t N_force_steps, uint8_t *status);

static void get_dFaux_increment(double *dFaux, double *dF_increment, double *dF_basic, uint32_t F_memsize, uint32_t *N_plastic_elem,
                        double *total_displacement, double *displacement, uint32_t N_nod);

static void get_plastic_elements(uint32_t N_elem, nb_plastified_analysis2D *elem_regime, bool *plastic_elements);

static int adjust_plastic_elem_to_yield_stress(double stress_tolerance, double *displacement, double *total_displacement,
                                        uint32_t F_memsize, uint32_t F_elemsize, uint32_t N_nod, uint32_t N_elem, double yield_stress,
                                        double *dF_increment, double *dFaux, double *max_vm_stress, uint32_t *plastified_elem,
                                        const nb_material_t *const material, nb_analysis2D_t analysis2D, nb_plastified_analysis2D *elem_regime,
                                        double *elastic_strain, const nb_fem_elem_t *const elemtype,
                                        const bool* elements_enabled, nb_sparse_t* K, const nb_mesh2D_t *const part,
                                        double *stress, uint32_t *N_plastic_elem, nb_analysis2D_params * params2D,
                                        uint32_t *simultaneous_elements, uint32_t N_simultaneous_plastic_elem, double *dF_basic,
                                        double *strain);

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
			 double *total_strain, /* Output */
			 double *stress, /*Output*/
			 double *total_displacement, /* Output */
			 uint32_t N_force_steps,
			 double stress_tolerance,
			 bool *plastic_elements)
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

    	uint32_t displacement_size = F_memsize;
    	uint32_t elastic_strain_size = F_elemsize;
    	uint32_t elem_regime_size = 10 * N_elem * sizeof(char);
    	uint32_t dF_size = F_memsize;
    	uint32_t dFaux_size = F_memsize;
    	uint32_t dF_basic_size = F_memsize;
    	uint32_t dF_increment_size = F_memsize;
    	uint32_t F_size = F_memsize;
    	uint32_t strain_size = F_elemsize;

    	uint64_t memsize = displacement_size + elastic_strain_size + elem_regime_size + dF_size + dFaux_size + dF_basic_size +
                       		 dF_increment_size + F_size + strain_size;

    	char* memblock = malloc(memsize);

    	double *displacement = (void*)memblock; //nb_allocate_mem(F_memsize); //
    	double *elastic_strain = (void*)(memblock + displacement_size); //nb_allocate_mem(F_elemsize); //
   	 nb_plastified_analysis2D *elem_regime = nb_allocate_mem(10 * N_elem * sizeof(char)); //(void*)(memblock + displacement_size + 			elastic_strain_size);
	double *dF = (void*)(memblock + displacement_size + elastic_strain_size + elem_regime_size); //nb_allocate_mem(F_memsize); //
    	double *dFaux = (void*)(memblock + displacement_size + elastic_strain_size + elem_regime_size + dF_size); //nb_allocate_mem(F_memsize); //
   	double *dF_basic = (void*)(memblock + displacement_size + elastic_strain_size + elem_regime_size + dF_size + dFaux_size); //nb_allocate_mem(F_memsize); //
    	double *dF_increment = (void*)(memblock + displacement_size + elastic_strain_size + elem_regime_size + dF_size + dFaux_size +
                                   dF_basic_size); //nb_allocate_mem(F_memsize);
   	double *F = (void*)(memblock + displacement_size + elastic_strain_size + elem_regime_size + dF_size + dFaux_size +
                        dF_basic_size + dF_increment_size); //nb_allocate_mem(F_memsize);
	double *strain = (void*)(memblock + displacement_size + elastic_strain_size + elem_regime_size + dF_size + dFaux_size +
                             dF_basic_size + dF_increment_size + F_size); //nb_allocate_mem(F_elemsize);
	memset(F, 0, F_memsize);
	memset(total_strain, 0, F_elemsize);
	memset(displacement, 0, F_memsize);
	memset(total_displacement, 0, F_memsize);

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

   	get_elem_regime(elem_regime, N_elem);

    	get_dF_basic(dF_basic, F_memsize, N_nod, F, N_force_steps);

    	memset(dF, 0, F_memsize);

    	double yield_stress = nb_material_get_yield_stress(material);
    	int iteration = 0;
    	uint8_t first_plastic_solver_status = 0;
    	iteration = get_first_plastic_element(K, dF_basic, displacement, N_nod, total_displacement, part, strain,
                                          elemtype, material, analysis2D, elements_enabled, stress, elem_regime,
                                          elastic_strain, iteration, N_force_steps, &first_plastic_solver_status);

   	if(0 != first_plastic_solver_status){
        	status = 2;
        	goto CLEANUP_LINEAR_SYSTEM;
   	}

    	for (int i = iteration; i < N_force_steps; i++) {
        	double max_vm_stress = 0;
        	uint32_t N_plastic_elem = 0;
        	uint32_t plastified_elem = 0;
        	uint32_t N_simultaneous_plastic_elem = 0;
        	printf("Iteration: %d\n", i); /* TEMPORAL */
        	int solver_status = plastic_solver(K, dF_basic, displacement, N_nod);

        	if(0 != solver_status){
            		status = 2;
            		goto CLEANUP_LINEAR_SYSTEM;
        	}

       		add_displacements(total_displacement, displacement, N_nod);
        	pipeline_compute_strain(strain, part, total_displacement, elemtype);
        	nb_fem_compute_plastic_stress_from_strain(N_elem, elemtype, material, analysis2D,
                                                      strain, elements_enabled, stress, elem_regime, elastic_strain);



        	get_stress_params(&max_vm_stress, &N_plastic_elem, &plastified_elem, &N_simultaneous_plastic_elem,
                          stress, yield_stress, elem_regime, N_elem, stress_tolerance);

        	uint32_t *simultaneous_elements = nb_allocate_mem(N_simultaneous_plastic_elem * sizeof(N_simultaneous_plastic_elem)); //= 			nb_soft_allocate_mem(N_simultaneous_plastic_elem*sizeof(uint32_t));

        	get_simultaneous_plastic_elements(simultaneous_elements, &max_vm_stress, &N_plastic_elem, N_simultaneous_plastic_elem,
                                          N_elem, stress_tolerance, stress);
       		printf("Simultaneous elements: %d\n", N_simultaneous_plastic_elem);

        	get_dFaux_increment(dFaux, dF_increment, dF_basic, F_memsize, &N_plastic_elem, total_displacement, displacement, N_nod);
	
        	uint8_t adjustment_to_yield_stress_status = adjust_plastic_elem_to_yield_stress(stress_tolerance, displacement, 											total_displacement, F_memsize, F_elemsize, N_nod,N_elem, 												yield_stress, dF_increment,dFaux,&max_vm_stress,
											&plastified_elem, material,analysis2D,elem_regime,
						                                     	elastic_strain, elemtype, elements_enabled, K, part, 												stress,&N_plastic_elem, params2D, simultaneous_elements,
											N_simultaneous_plastic_elem, dF_basic, strain);
       		if (adjustment_to_yield_stress_status != 0)
            		goto CLEANUP_LINEAR_SYSTEM;
    	}

    	add_total_strain_plastic(strain, total_strain, N_elem);

    	get_plastic_elements(N_elem, elem_regime, plastic_elements);

    	//print_results_on_graph(N_nod, N_elem, total_displacement, stress, elem_regime, total_strain, part, plastic_elements);

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
			 nb_plastified_analysis2D *elem_regime,
			 double *elastic_strain)
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
                        			stress[id * 3] = (strain[id * 3] - elastic_strain[id * 3]) * Dp[0] +
                            			(strain[id*3 + 1] - elastic_strain[id*3 + 1]) * Dp[1] +
                            			elastic_strain[id * 3]*De[0] + elastic_strain[id*3 + 1]*De[1];
                        			stress[id*3 + 1] = (strain[id * 3] - elastic_strain[id * 3]) * Dp[1] +
                            			(strain[id*3 + 1] - elastic_strain[id*3 + 1]) * Dp[2] +
                            			elastic_strain[id*3 + 1]*De[1] + elastic_strain[id*3 + 1]*De[2];
                       				stress[id*3 + 2] = (strain[id*3 + 2] - elastic_strain[id*3 + 2]) * Dp[3] +
                       				elastic_strain[id*3 + 2]*De[3];
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

void add_displacements(double *total_displacement, double *displacement, uint32_t N_nod) {
   	for(int i = 0; i < 2*N_nod; i++) {
        	total_displacement[i] += displacement[i];
    	}
}

void substract_displacement(double *total_displacement, double *displacement, uint32_t N_nod) {
    	for(int i = 0; i < 2*N_nod; i++) {
        	total_displacement[i] -= displacement[i];
    	}
}

int add_elastic_displacements(double *total_displacement, double* displacement, uint32_t N_nod,
                              const nb_mesh2D_t *const part,double *strain, const nb_fem_elem_t *const elemtype,
                              const nb_material_t *const material, nb_analysis2D_t analysis2D,
                              const bool* elements_enabled,double *stress, nb_plastified_analysis2D *elem_regime,
                              double *elastic_strain, int i, uint32_t N_force_steps)
{
    	uint32_t N_elem = nb_mesh2D_get_N_elems(part);
    	double yield_stress = nb_material_get_yield_stress(material);
    	double VM_stress = 0;
    	while(i < N_force_steps) {
        	add_displacements(total_displacement, displacement, N_nod);
        	pipeline_compute_strain(strain, part, total_displacement, elemtype);
        	nb_fem_compute_plastic_stress_from_strain(N_elem, elemtype, material, analysis2D,
                                                    strain, elements_enabled, stress, elem_regime, elastic_strain);
        	for(int j = 0; j < N_elem; j++){
            		VM_stress = nb_pde_get_vm_stress(stress[3*j], stress[3*j +1], stress[3*j +2]);
            		if(VM_stress > yield_stress){
                		printf("VM[%d]: %lf\n", j, VM_stress); /* TEMPORAL */
                		substract_displacement(total_displacement, displacement, N_nod);
                		return i;
            		}
        	}
        	printf("Iteration: %d\n", i); /* TEMPORAL */
        	i++;
    	}
    	return i;
}

void get_stress_params(double *max_vm_stress, uint32_t *N_plastic_elem,
                       uint32_t *plastified_elem, uint32_t *N_simultaneous_plastic_elem, double *stress, double yield_stress,
                       nb_plastified_analysis2D *elem_regime,
                       uint32_t N_elem, double stress_tolerance) {
    	double vm_stress = 0;
    	max_vm_stress[0] = 0;    //memset(max_vm_stress, 0, sizeof(double));
    	N_plastic_elem[0] = 0;
    	//memset(N_plastic_elem, 0, sizeof(uint32_t));
    	N_plastic_elem[0] = 1;
    	for(int j = 0; j < N_elem; j++) {
        	vm_stress = nb_pde_get_vm_stress(stress[3*j], stress[3*j +1], stress[3*j +2]);
        	if(vm_stress > yield_stress && elem_regime[j] != NB_PLASTIC) {
            		N_plastic_elem[0] += 1;
            		if(vm_stress > max_vm_stress[0]) {
                		max_vm_stress[0] = vm_stress;
                		plastified_elem[0] = j;
            		}
        	}
    	}
   	N_simultaneous_plastic_elem[0] = 0;
    	//memset(N_simultaneous_plastic_elem, 0, sizeof(uint32_t));
    	for(int j = 0; j < N_elem; j++) {
        	vm_stress = nb_pde_get_vm_stress(stress[3*j], stress[3*j + 1], stress[3*j + 2]);
        	if((vm_stress - max_vm_stress[0]) < 0 && (vm_stress - max_vm_stress[0]) > -stress_tolerance){
           		N_simultaneous_plastic_elem[0] += 1;
            		N_plastic_elem[0] -=1;
        	}
        	else if((vm_stress - max_vm_stress[0]) >= 0 && (vm_stress - max_vm_stress[0]) < stress_tolerance){
            		N_simultaneous_plastic_elem[0] += 1;
           		N_plastic_elem[0] -=1;
        	}
    	}

    	if(N_simultaneous_plastic_elem[0] == 0){
        	N_plastic_elem[0] = 0;
    	}
}

void get_simultaneous_plastic_elements(uint32_t *simultaneous_elements,
                                        double *max_vm_stress, uint32_t *N_plastic_elem,
                                        uint32_t N_simultaneous_plastic_elem,
                                        uint32_t N_elem,
                                        double stress_tolerance,
                                        double *stress) 
{
    	uint32_t sim_elem = -1;
    	double vm_stress;
    	if(N_simultaneous_plastic_elem > 1) {
        	for(int k = 0; k < N_elem; k++) {
            		vm_stress = nb_pde_get_vm_stress(stress[3*k], stress[3*k +1], stress[3*k +2]);
            		if((vm_stress - max_vm_stress[0]) < 0 && (vm_stress - max_vm_stress[0]) > -stress_tolerance){
                		sim_elem += 1;
               			simultaneous_elements[sim_elem] = k;
            		}
            		else if((vm_stress - max_vm_stress[0]) >= 0 && (vm_stress - max_vm_stress[0]) < stress_tolerance){
                		sim_elem += 1;
                		simultaneous_elements[sim_elem] = k;
            		}
        	}
    	}
 }

void find_new_plastic_elem(uint32_t N_elem,
                           nb_plastified_analysis2D *elem_regime,
                           double *stress, double yield_stress,
                           double *max_vm_stress,
                           uint32_t *plastified_elem,
                           uint32_t *N_plastic_elem)
{
    	uint32_t N_plastic_elem_aux = 0;
    	plastified_elem[0] = 0;
    	max_vm_stress[0] = 0;
    	//memset(plastified_elem, 0, sizeof(uint32_t));
    	//memset(max_vm_stress, 0, sizeof(double));
    	double vm_stress;

    	for(int j = 0; j < N_elem; j++) {
        	vm_stress = nb_pde_get_vm_stress(stress[3*j], stress[3*j +1], stress[3*j +2]);
        	if(vm_stress > yield_stress && elem_regime[j] != NB_PLASTIC) {
            		N_plastic_elem_aux += 1;
            		if(vm_stress > max_vm_stress[0] && elem_regime[j] != NB_PLASTIC) {
                		max_vm_stress[0] = vm_stress;
                		plastified_elem[0] = j;
            		}
        	}
    	}

    	if(N_plastic_elem_aux > N_plastic_elem[0]){
        	N_plastic_elem[0] = N_plastic_elem_aux;
    	}
}

int adjust_force_to_yield_stress(double limit_stress, double stress_tolerance, double *displacement, double *total_displacement,
                                     uint32_t F_memsize, uint32_t F_elemsize, uint32_t N_nod, uint32_t N_elem, double yield_stress,
                                     double *dF_increment, double *dFaux, double *max_vm_stress, uint32_t *plastified_elem,
                                     const nb_material_t *const material, nb_analysis2D_t analysis2D, nb_plastified_analysis2D *elem_regime,
                                     double *elastic_strain, const nb_fem_elem_t *const elemtype, double *aux_strain,
                                     const bool* elements_enabled, const nb_sparse_t *const K, const nb_mesh2D_t *const part,
                                     double *stress)
{
    	uint32_t iter = 0;
    	while(limit_stress < -stress_tolerance || limit_stress > stress_tolerance){
        	memset(displacement, 0, F_memsize);
        	memset(aux_strain, 0, F_elemsize);
        	if(iter > 0 && max_vm_stress[0] > yield_stress) {
            		for(int k = 0; k < 2*N_nod; k++) {
                		if(iter == 1)
                    			dF_increment[k] = dFaux[k]/2;
                		else
                    			dF_increment[k] /= 2;

                		dFaux[k] -= dF_increment[k];
            		}
        	}

        	if(iter > 0 && max_vm_stress[0] < yield_stress) {
            		for(int k = 0; k < 2*N_nod; k++) {
                		if(iter == 1)
                    			dF_increment[k] = dFaux[k]/2;
                		else
                    			dF_increment[k] /= 2;
               			dFaux[k] += dF_increment[k];
            		}
        	}

        	iter += 1;

        	int solver_status = plastic_solver(K, dFaux, displacement, N_nod);
        	if(0 != solver_status)
            		return 2;
        	
        	add_displacements(total_displacement, displacement, N_nod);

        	pipeline_compute_strain(aux_strain, part, total_displacement, elemtype);

        	substract_displacement(total_displacement, displacement, N_nod);

       		nb_fem_compute_plastic_stress_from_strain(N_elem, elemtype, material, analysis2D,
                                                              aux_strain, elements_enabled, stress, elem_regime, elastic_strain);

        	max_vm_stress[0] = nb_pde_get_vm_stress(stress[3*plastified_elem[0]], 
							stress[3*plastified_elem[0] + 1], 
							stress[3*plastified_elem[0] + 2]);

       		limit_stress = max_vm_stress[0] - yield_stress;
    	}
    	return 0;
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

void add_total_strain_plastic(double *strain, double *total_strain, uint32_t N_elem) {
    	for (int i = 0; i < N_elem; i++) {
        	total_strain[3*i] = strain[3*i];
        	total_strain[3*i + 1] = strain[3*i + 1];
        	total_strain[3*i + 2] = strain[3*i + 2];
    	}
}

static void get_elem_regime(nb_plastified_analysis2D *elem_regime, uint32_t N_elem){
    	for (int i = 0; i < N_elem; i++) {
        	elem_regime[i] = NB_ELASTIC;
    	}
}

static void get_dF_basic(double *dF_basic, uint32_t F_memsize, uint32_t N_nod, double *F, uint32_t N_force_steps) {
    	memset(dF_basic, 0, F_memsize);
	for (int i = 0; i < 2*N_nod; i++) {
        	dF_basic[i] = F[i] / N_force_steps;
    	}
}

static int get_first_plastic_element(const nb_sparse_t *const K, double *dF_basic, double *displacement, uint32_t N_nod,
                              double *total_displacement, const nb_mesh2D_t *const part, double *strain,
                              const nb_fem_elem_t *const elemtype, const nb_material_t *const material,
                              nb_analysis2D_t analysis2D, const bool *elements_enabled, double *stress,
                              nb_plastified_analysis2D *elem_regime, double *elastic_strain,
                              uint32_t i, uint32_t N_force_steps, uint8_t *status)
{
    	status[0] = plastic_solver(K, dF_basic, displacement, N_nod);

    	i = add_elastic_displacements(total_displacement, displacement, N_nod, part, strain, elemtype,
                                	material, analysis2D, elements_enabled, stress, elem_regime,
                                	elastic_strain, i, N_force_steps);
   	 return i;
}

static void get_dFaux_increment(double *dFaux, double *dF_increment, double *dF_basic, uint32_t F_memsize, uint32_t *N_plastic_elem,
                    double *total_displacement, double *displacement,uint32_t N_nod)
{
    	memset(dFaux, 0, F_memsize);
    	memset(dF_increment, 0, F_memsize);
    	if(N_plastic_elem[0] != 0) {
       		substract_displacement(total_displacement, displacement, N_nod);
        	for(int k = 0; k < 2*N_nod; k++){
            		dFaux[k] = dF_basic[k]/2;
        	}
   	}
}

static int adjust_plastic_elem_to_yield_stress(double stress_tolerance, double *displacement, double *total_displacement,
                                        uint32_t F_memsize, uint32_t F_elemsize, uint32_t N_nod, uint32_t N_elem, double yield_stress,
                                        double *dF_increment, double *dFaux, double *max_vm_stress, uint32_t *plastified_elem,
                                        const nb_material_t *const material, nb_analysis2D_t analysis2D, nb_plastified_analysis2D *elem_regime,
                                        double *elastic_strain, const nb_fem_elem_t *const elemtype,
                                        const bool* elements_enabled, nb_sparse_t* K, const nb_mesh2D_t *const part,
                                        double *stress, uint32_t *N_plastic_elem, nb_analysis2D_params * params2D,
                                        uint32_t *simultaneous_elements, uint32_t N_simultaneous_plastic_elem, double *dF_basic,
                                        double *strain)
{
    	int status = 0;
    	while(N_plastic_elem[0] != 0) {
        	double limit_stress = 0;
        	max_vm_stress[0] = 0;
        	//memset(max_vm_stress, 0, sizeof(double));
        	max_vm_stress[0] = nb_pde_get_vm_stress(stress[3*plastified_elem[0]], stress[3*plastified_elem[0] + 1], 
							stress[3*plastified_elem[0] + 2]);
        	limit_stress = max_vm_stress[0] - yield_stress;
        	double *aux_strain = nb_allocate_mem(F_elemsize);
        	printf("Plastified elements: %d\n", N_plastic_elem[0]);

        	int adjustment_status = adjust_force_to_yield_stress(limit_stress, stress_tolerance, displacement, total_displacement,
                                                                    F_memsize, F_elemsize, N_nod, N_elem, yield_stress, dF_increment,
                                                                    dFaux, max_vm_stress, plastified_elem, material, analysis2D, elem_regime,
                                                                    elastic_strain, elemtype, aux_strain, elements_enabled, K, part, stress);
        	if(adjustment_status != 0){
            		status = 3;
            		goto EXIT;
        	}
        	add_displacements(total_displacement, displacement, N_nod);

        	int update_status = updated_stiffness_matrix(K, elem_regime,
                                                    elements_enabled, part,
                                                    elemtype, material,
                                                    analysis2D, params2D, elastic_strain,
                                                    simultaneous_elements, aux_strain,
                                                    N_simultaneous_plastic_elem,
                                                    plastified_elem,
                                                    N_plastic_elem);
        	if(update_status != 0){
            		status = 2;
            		goto EXIT;
        	}

        	for (int k = 0; k < 2*N_nod; k++) {
            		dFaux[k] = dF_basic[k] - dFaux[k];
        	}

        	int solver_status = plastic_solver(K, dFaux, displacement, N_nod);

        	if(0 != solver_status){
            		status = 2;
            		goto EXIT;
        	}
        	add_displacements(total_displacement, displacement, N_nod);

        	pipeline_compute_strain(strain, part, total_displacement, elemtype);

        	nb_fem_compute_plastic_stress_from_strain(N_elem, elemtype, material, analysis2D,
                                                    strain, elements_enabled, stress, elem_regime, elastic_strain);

        	find_new_plastic_elem(N_elem, elem_regime, stress, yield_stress, max_vm_stress, plastified_elem, N_plastic_elem);

        	if (N_plastic_elem[0] != 0) {
            		substract_displacement(total_displacement, displacement, N_nod);
            		for (int k = 0; k < 2*N_nod; k++) {
                		dFaux[k] = dFaux[k]/2;
            		}
        	}

    		nb_free_mem(aux_strain);
    	}
    	EXIT:
    	return status;
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
