#ifndef STATIC_PLASTICITY2D_H_INCLUDED
#define STATIC_PLASTICITY2D_H_INCLUDED

#include "nb/solver_bot.h"
#include "nb/geometric_bot.h"
#include "nb/pde_bot/material.h"
#include "nb/pde_bot/boundary_conditions/bcond.h"
#include "nb/pde_bot/finite_element/element.h"
#include "../../../../../sources/nb/pde_bot/finite_element/solid_mechanics/pipeline.h"
#include "nb/pde_bot/common_solid_mechanics/analysis2D.h"
#include "nb/pde_bot/common_solid_mechanics/formulas.h"

void add_displacements(double *total_displacement, double *displacement, uint32_t N_nod);
void substract_displacement(double *total_displacement, double *displacement, uint32_t N_nod);
void add_elastic_strain(double *strain, double *elastic_strain, nb_plastified_analysis2D *elem_regime, uint32_t N_elem);
void add_plastic_strain(double *strain, double *plastic_strain, nb_plastified_analysis2D *elem_regime, uint32_t N_elem);
void add_total_strain(double *total_strain, double *elastic_strain, double *plastic_strain, uint32_t N_elem);
void add_total_strain_plastic(double *strain, double *total_strain, uint32_t N_elem);
void add_internal_forces_to_global_system(const nb_fem_elem_t *elem, uint32_t id,
                const nb_mesh2D_t *part, double *FI,
                double *FIe);
int assemble_internal_forces_element(const nb_fem_elem_t *elem,
                uint32_t id, double *stress,
			    const nb_mesh2D_t *part,
			    const nb_material_t *material,
			    nb_analysis2D_t analysis2D,
			    nb_analysis2D_params *params2D,
			    uint32_t N_nod,
			    double *FI);
int compute_internal_forces
                (double *FI,
                double *stress,
                uint32_t N_elem,
                uint32_t N_nod,
                const nb_mesh2D_t *part,
                const nb_material_t *material,
                nb_analysis2D_t analysis2D,
                nb_analysis2D_params *params2D,
                const nb_fem_elem_t *const elem);
int integrate_elemental_FIe_vector
                (const nb_fem_elem_t *elem, uint32_t id,
                const nb_mesh2D_t *part, double *stress,
                nb_analysis2D_params *params2D,
                uint32_t N_nod,
                double *FIe);
void internal_forces_sum_gauss_point(const nb_fem_elem_t *elem, uint32_t gp_id, uint8_t GP,
                double thickness, double detJ, const nb_mesh2D_t *part,
                double *dNi_dx, double *dNi_dy,
                double *stress, double *FIe);
double force_tolerance (double *F, double *FI, uint32_t N_nod);
void nb_fem_compute_elastic_stress_from_strain
			(uint32_t N_elements,
			 const nb_fem_elem_t *const elem,
			 const nb_material_t *const material,
			 nb_analysis2D_t analysis2D,
			 double* strain,
			 const bool* elements_enabled,
			 double* stress);
void nb_fem_compute_plastic_stress_from_strain_IF
                (uint32_t N_elements,
                const nb_fem_elem_t *const elem,
                const nb_material_t *const material,
                nb_analysis2D_t analysis2D,
                double* strain,
                const bool* elements_enabled /* NULL to enable all */,
                double* stress /* Output */,
                nb_plastified_analysis2D *elem_regime);
void nb_fem_compute_plastic_stress_from_strain
			(uint32_t N_elements,
			 const nb_fem_elem_t *const elem,
			 const nb_material_t *const material,
			 nb_analysis2D_t analysis2D,
			 double* strain,
			 const bool* elements_enabled ,
			 double* stress ,
			 nb_plastified_analysis2D *elem_regime,
			 double *elastic_strain);
void nb_fem_compute_diference_of_stresses(double *residual_stress, double *elastic_stress,
                                          double *plastic_stress, uint32_t N_elements,
                                          nb_plastified_analysis2D *elem_regime);
void residual_force(double *res_force, double *dF, double *FI, uint32_t N_nod);
void nb_compute_plastic_strain(double *plastic_strain, double *total_strain, double *last_elastic_strain,
                               uint32_t N_elem, nb_plastified_analysis2D *elem_regime);
int plastic_solver(const nb_sparse_t *const A, const double *const b, double* x, uint32_t N_nod);
double stress_max_tolerance(double *plastic_stress, nb_plastified_analysis2D *elem_regime, uint32_t N_elements);
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
			 double accepted_tol,
			 bool *plastic_elements);
double stress_max_tolerance(double *plastic_stress, nb_plastified_analysis2D *elem_regime, uint32_t N_elem);
int add_elastic_displacements(double *total_displacement, double* displacement, uint32_t N_nod, const nb_mesh2D_t *const part,
                              double *strain, const nb_fem_elem_t *const elemtype, const nb_material_t *const material,
                              nb_analysis2D_t analysis2D, const bool* elements_enabled,double *stress,
                              nb_plastified_analysis2D *elem_regime, double *elastic_strain,
                              int i, uint32_t N_force_steps);
void get_stress_params(double *max_vm_stress, uint32_t *N_plastic_elem,
                       uint32_t *plastified_elem, uint32_t *N_simultaneous_plastic_elem, double *stress, double yield_stress,
                       nb_plastified_analysis2D *elem_regime,
                       uint32_t N_elem, double stress_tolerance);
 void get_simultaneous_plastic_elements(uint32_t *simultaneous_elements,
                                        double *max_vm_stress, uint32_t *N_plastic_elem,
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
int adjust_force_to_yield_stress(double limit_stress, double stress_tolerance, double *displacement, double *total_displacement,
                                     uint32_t F_memsize, uint32_t F_elemsize, uint32_t N_nod, uint32_t N_elem, double yield_stress,
                                     double *dF_increment, double *dFaux, double *max_vm_stress, uint32_t *plastified_elem,
                                     const nb_material_t *const material, nb_analysis2D_t analysis2D, nb_plastified_analysis2D *elem_regime,
                                     double *elastic_strain, const nb_fem_elem_t *const elemtype, double *aux_strain,
                                     const bool* elements_enabled, const nb_sparse_t *const K, const nb_mesh2D_t *const part,
                                     double *stress);

#endif // STATIC_PLASTICITY2D_H_INCLUDED
