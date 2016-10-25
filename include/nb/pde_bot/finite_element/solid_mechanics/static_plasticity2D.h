#ifndef STATIC_PLASTICITY2D_H_INCLUDED
#define STATIC_PLASTICITY2D_H_INCLUDED

#include "nb/solver_bot.h"
#include "nb/geometric_bot.h"
#include "nb/pde_bot/material.h"
#include "nb/pde_bot/boundary_conditions/bcond.h"
#include "nb/pde_bot/finite_element/element.h"
#include "nb/pde_bot/finite_element/solid_mechanics/pipeline.h"
#include "nb/pde_bot/common_solid_mechanics/analysis2D.h"
#include "nb/pde_bot/common_solid_mechanics/formulas.h"

void add_elastic_strain(double *strain, double *elastic_strain, nb_plastified_analysis2D *elem_regime, uint32_t N_elem);
void add_plastic_strain(double *strain, double *plastic_strain, uint32_t N_elem);
void add_total_strain(double *total_strain, double *elastic_strain, double *plastic_strain, uint32_t N_elem);
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
void internal_forces_sum_gauss_point(const nb_fem_elem_t *elem, uint32_t gp_id,
                double thickness, double detJ, const nb_mesh2D_t *part,
                double *dNi_dx, double *dNi_dy,
                double *stress, double *FIe);
double force_tolerance (double *F, double *FI, uint32_t N_nod);
void nb_fem_compute_plastic_stress_from_strain
                (uint32_t N_elements,
                const nb_fem_elem_t *const elem,
                const nb_material_t *const material,
                nb_analysis2D_t analysis2D,
                double* strain,
                const bool* elements_enabled /* NULL to enable all */,
                double* stress /* Output */,
                nb_plastified_analysis2D *elem_regime);
void residual_force(double *res_force, double *dF, double *FI, uint32_t N_nod);
int plastic_solver(const nb_sparse_t *const A, const double *const b, double* x);
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
			 double N_force_steps,
			 double accepted_tol);

#endif // STATIC_PLASTICITY2D_H_INCLUDED
