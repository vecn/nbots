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
int add_elastic_displacements(double *total_displacement, double* displacement, uint32_t N_nod, const nb_mesh2D_t *const part,
                              double *strain, const nb_fem_elem_t *const elemtype, const nb_material_t *const material,
                              nb_analysis2D_t analysis2D, const bool* elements_enabled,double *stress,
                              nb_plastified_analysis2D *elem_regime, double *elastic_strain,
                              bool elastic_regime, int i, uint32_t N_force_steps);
void get_stress_params(double *max_vm_stress, uint32_t *N_plastic_elem,
                       uint32_t *plastified_elem, uint32_t *N_simultaneous_plastic_elem, double *stress, double yield_stress,
                       nb_plastified_analysis2D *elem_regime,
                       uint32_t N_elem, double stress_tolerance);
 void get_simultaneous_plastic_elements(uint32_t *simultaneous_elements,
                                        double *max_vm_stress, uint32_t N_plastic_elem,
                                        uint32_t *N_simultaneous_plastic_elem,
                                        uint32_t N_elem,
                                        double stress_tolerance,
                                        double *stress);
void find_new_plastic_elem(uint32_t N_elem,
                           nb_plastified_analysis2D *elem_regime,
                           double *stress, double yield_stress,
                           double *max_vm_stress,
                           uint32_t *plastified_elem,
                           uint32_t *N_plastic_elem);
uint8_t adjust_force_to_yield_stress(double limit_stress, double stress_tolerance, double *displacement, double *total_displacement,
                                     uint32_t F_memsize, uint32_t F_elemsize, uint32_t N_nod, uint32_t N_elem, double yield_stress,
                                     double *dF_increment, double *dFaux, double *max_vm_stress, uint32_t *plastified_elem,
                                     const nb_material_t *const material, nb_analysis2D_t analysis2D, nb_plastified_analysis2D *elem_regime,
                                     double *elastic_strain, const nb_fem_elem_t *const elemtype, double *aux_strain,
                                     const bool* elements_enabled, const nb_sparse_t *const K, const nb_mesh2D_t *const part,
                                     double *stress);

int fem_compute_plastic_2D_Solid_Mechanics_bis
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
			 double *stress,
			 double *total_displacement, /* Output */
			 uint32_t N_force_steps,
			 double stress_tolerance)
{
	int status = 0;
	nb_graph_t *graph = malloc(nb_graph_get_memsize());
	nb_graph_init(graph);
	nb_mesh2D_load_graph(part, graph, NB_NODES_LINKED_BY_ELEMS);
	nb_sparse_t *K = nb_sparse_create(graph, NULL, 2);
	nb_graph_finish(graph);

	uint32_t N_nod = nb_mesh2D_get_N_nodes(part);
	uint32_t N_elem = nb_mesh2D_get_N_elems(part);
	uint32_t F_memsize = 2 * N_nod * sizeof(double);
	uint32_t F_elemsize = 3*N_elem*sizeof(double);
	double *F = nb_soft_allocate_mem(F_memsize);
	memset(F, 0, F_memsize);
	double *strain = nb_allocate_mem(F_elemsize);
	double *displacement = nb_allocate_mem(F_memsize);
	double *elastic_strain = nb_allocate_mem(F_elemsize);
	memset(total_strain, 0, F_elemsize);

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

    double *dF_basic = nb_soft_allocate_mem(F_memsize);
    memset(dF_basic, 0, F_memsize);
	for (int i = 0; i < 2*N_nod; i++) {
        dF_basic[i] = F[i] / N_force_steps;
    }

    double *dF = nb_soft_allocate_mem(F_memsize);
    memset(dF, 0, F_memsize);
    double *dFaux = nb_soft_allocate_mem(F_memsize);
    double *dF_increment = nb_soft_allocate_mem(F_memsize);
    bool elastic_regime = true;
    double yield_stress = nb_material_get_yield_stress(material);

    for (int i = 0; i < N_force_steps; i++) {

        if(elastic_regime){
            int solver_status = plastic_solver(K, dF_basic, displacement);
            if(0 != solver_status){
                status = 2;
                goto CLEANUP_LINEAR_SYSTEM;
            }
            i = add_elastic_displacements(total_displacement, displacement, N_nod, part, strain, elemtype,
                                          material, analysis2D, elements_enabled, stress, elem_regime,
                                          elastic_strain, elastic_regime, i, N_force_steps);
            elastic_regime = false;
        }
        printf("Iteration: %d\n", i); /* TEMPORAL */
        int solver_status = plastic_solver(K, dF_basic, displacement);
        if(0 != solver_status){
            status = 2;
            goto CLEANUP_LINEAR_SYSTEM;
        }

        add_displacements(total_displacement, displacement, N_nod);
        pipeline_compute_strain(strain, part, total_displacement, elemtype);
        nb_fem_compute_plastic_stress_from_strain_bis(N_elem, elemtype, material, analysis2D,
                                                      strain, elements_enabled, stress, elem_regime, elastic_strain);

        double *max_vm_stress = nb_allocate_mem(sizeof(double));
        uint32_t *N_plastic_elem = nb_allocate_mem(sizeof(uint32_t));
        uint32_t *plastified_elem = nb_allocate_mem(sizeof(uint32_t));
        uint32_t *N_simultaneous_plastic_elem = nb_allocate_mem(sizeof(uint32_t));

        get_stress_params(max_vm_stress, N_plastic_elem, plastified_elem, N_simultaneous_plastic_elem,
                          stress, yield_stress, elem_regime, N_elem, stress_tolerance);

        uint32_t *simultaneous_elements = nb_soft_allocate_mem(N_simultaneous_plastic_elem[0]*sizeof(uint32_t));

        get_simultaneous_plastic_elements(simultaneous_elements, max_vm_stress, N_plastic_elem, N_simultaneous_plastic_elem,
                                          N_elem, stress_tolerance, stress);

        memset(dFaux, 0, F_memsize);
        memset(dF_increment, 0, F_memsize);
        if(N_plastic_elem[0] != 0) {
            substract_displacement(total_displacement, displacement, N_nod);
            for(int k = 0; k < 2*N_nod; k++){
                dFaux[k] = dF_basic[k]/2;
            }
        }

        while(N_plastic_elem[0] != 0) {
            double limit_stress = 0;
            memset(max_vm_stress, 0, sizeof(double));
            max_vm_stress[0] = nb_pde_get_vm_stress(stress[3*plastified_elem[0]], stress[3*plastified_elem[0] + 1], stress[3*plastified_elem[0] + 2]);
            limit_stress = max_vm_stress[0] - yield_stress;

            double *aux_strain = nb_allocate_mem(F_elemsize);
            printf("Plastified elements: %d\n", N_plastic_elem[0]);

            uint32_t adjustment_status = adjust_force_to_yield_stress(limit_stress, stress_tolerance, displacement, total_displacement,
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
               goto EXIT;
            }

            for (int k = 0; k < 2*N_nod; k++) {
                dFaux[k] = dF_basic[k] - dFaux[k];
            }

            int solver_status = plastic_solver(K, dFaux, displacement);

            if(0 != solver_status){
                status = 2;
                goto CLEANUP_LINEAR_SYSTEM;
            }
            add_displacements(total_displacement, displacement, N_nod);

            pipeline_compute_strain(strain, part, total_displacement, elemtype);

            nb_fem_compute_plastic_stress_from_strain_bis(N_elem, elemtype, material, analysis2D,
                                                      strain, elements_enabled, stress, elem_regime, elastic_strain);

            find_new_plastic_elem(N_elem, elem_regime, stress, yield_stress, max_vm_stress, plastified_elem, N_plastic_elem);

            if (N_plastic_elem[0] != 0) {
                substract_displacement(total_displacement, displacement, N_nod);
                for (int k = 0; k < 2*N_nod; k++) {
                    dFaux[k] = dFaux[k]/2;
                }
            }

            EXIT:
            nb_free_mem(aux_strain);
            if(status == 2 || update_status == 4 || status == 3){
                goto CLEANUP_LINEAR_SYSTEM;
            }
        }
        nb_free_mem(max_vm_stress);
        nb_free_mem(N_plastic_elem);
        nb_free_mem(plastified_elem);
        nb_free_mem(N_simultaneous_plastic_elem);
    }

    double *avg_displacement = nb_allocate_mem(N_nod*sizeof(double));
    for(int i = 1; i < N_nod; i++){
        avg_displacement[i] = sqrt(POW2(total_displacement[2*i])+ POW2(total_displacement[2*i +1]));
    }
    double *total_displacement_y = nb_allocate_mem(N_nod*sizeof(double));
    double *total_displacement_x = nb_allocate_mem(N_nod*sizeof(double));
    for(int j = 0; j < N_nod; j++){
        total_displacement_x[j] = total_displacement[2*j];
        total_displacement_y[j] = total_displacement[2*j +1];
    }
    double *Sx = nb_allocate_mem(N_elem*sizeof(double));
    double *Sy = nb_allocate_mem(N_elem*sizeof(double));
    double *Sxy = nb_allocate_mem(N_elem*sizeof(double));
    for(int j = 0; j < N_elem; j++) {
        Sx[j] = stress[3*j];
        Sy[j] = stress[3*j + 1];
        Sxy[j] = stress[3*j +2];
    }
    bool *plastic_elements = nb_allocate_mem(N_elem*sizeof(bool));
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
    nb_mesh2D_distort_with_field(part, NB_NODE, total_displacement, 0.5);
    nb_mesh2D_export_draw(part, "continuous_beam_plastic_elements.png", 1200, 800, NB_ELEMENT, NB_CLASS, plastic_elements, true);
    nb_mesh2D_export_draw(part,"continuous_beam_avg_disp.png", 1200, 800, NB_NODE, NB_FIELD, avg_displacement, true);
    nb_mesh2D_export_draw(part,"continuous_beam_disp_x.png", 1200, 800, NB_NODE, NB_FIELD, total_displacement_x, true);
    nb_mesh2D_export_draw(part,"continuous_beam_disp_y.png", 1200, 800, NB_NODE, NB_FIELD, total_displacement_y, true);
    nb_mesh2D_export_draw(part,"continuous_beam_Sx.png", 1200, 800, NB_ELEMENT, NB_FIELD, Sx, true);
    nb_mesh2D_export_draw(part,"continuous_beam_Sy.png", 1200, 800, NB_ELEMENT, NB_FIELD, Sy, true);
    nb_mesh2D_export_draw(part,"continuous_beam_Sxy.png", 1200, 800, NB_ELEMENT, NB_FIELD, Sxy, true);

    CLEANUP_LINEAR_SYSTEM:
    nb_free_mem(total_displacement_x);
    nb_free_mem(total_displacement_y);
    nb_free_mem(avg_displacement);
    nb_free_mem(Sx);
    nb_free_mem(Sy);
    nb_free_mem(Sxy);
    nb_sparse_destroy(K);
	nb_soft_free_mem(F_memsize, F);
	nb_soft_free_mem(F_memsize, dF_basic);
	nb_soft_free_mem(F_memsize, dF);
	nb_soft_free_mem(F_memsize, dFaux);
	nb_soft_free_mem(F_memsize, dF_increment);
	nb_free_mem(displacement);
	nb_free_mem(strain);
	nb_free_mem(elem_regime);

	return status;
}

void nb_fem_compute_plastic_stress_from_strain_bis
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

	//uint32_t omp_parallel_threads = 1;
    #pragma omp parallel for num_threads(omp_parallel_threads) schedule(guided)
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
                              double *elastic_strain, bool elastic_regime, int i, uint32_t N_force_steps) {
    uint32_t N_elem = nb_mesh2D_get_N_elems(part);
    double yield_stress = nb_material_get_yield_stress(material);
    double vm_stress;
    while(i < N_force_steps) {
        add_displacements(total_displacement, displacement, N_nod);
        pipeline_compute_strain(strain, part, total_displacement, elemtype);
        nb_fem_compute_plastic_stress_from_strain_bis(N_elem, elemtype, material, analysis2D,
                                                    strain, elements_enabled, stress, elem_regime, elastic_strain);
        for(int j = 0; j < N_elem; j++){
            vm_stress = nb_pde_get_vm_stress(stress[3*j], stress[3*j +1], stress[3*j +2]);
            if(vm_stress > yield_stress){
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
    double vm_stress;
    memset(max_vm_stress, 0, sizeof(double));
    memset(N_plastic_elem, 0, sizeof(uint32_t));
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
    memset(N_simultaneous_plastic_elem, 0, sizeof(uint32_t));
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
                                        double *max_vm_stress, uint32_t N_plastic_elem,
                                        uint32_t *N_simultaneous_plastic_elem,
                                        uint32_t N_elem,
                                        double stress_tolerance,
                                        double *stress) {
    uint32_t sim_elem = -1;
    double vm_stress;
    if(N_simultaneous_plastic_elem[0] > 1) {
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
    memset(plastified_elem, 0, sizeof(uint32_t));
    memset(max_vm_stress, 0, sizeof(double));
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

uint8_t adjust_force_to_yield_stress(double limit_stress, double stress_tolerance, double *displacement, double *total_displacement,
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

        int solver_status = plastic_solver(K, dFaux, displacement);
        if(0 != solver_status){
            return 2;
        }

        add_displacements(total_displacement, displacement, N_nod);

        pipeline_compute_strain(aux_strain, part, total_displacement, elemtype);

        substract_displacement(total_displacement, displacement, N_nod);

        nb_fem_compute_plastic_stress_from_strain_bis(N_elem, elemtype, material, analysis2D,
                                                              aux_strain, elements_enabled, stress, elem_regime, elastic_strain);

        max_vm_stress[0] = nb_pde_get_vm_stress(stress[3*plastified_elem[0]], stress[3*plastified_elem[0] + 1], stress[3*plastified_elem[0] + 2]);

        limit_stress = max_vm_stress[0] - yield_stress;
    }
    return 0;
}
