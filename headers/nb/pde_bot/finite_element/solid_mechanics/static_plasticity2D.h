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

void nb_fem_compute_plastic_stress_from_strain
			(uint32_t N_elements,
			 const nb_fem_elem_t *const elem,
			 const nb_material_t *const material,
			 nb_analysis2D_t analysis2D,
			 double* strain,
			 const bool* elements_enabled ,
			 double* stress ,
			 nb_plastified_analysis2D *elem_regime);

int plastic_solver(const nb_sparse_t *const A, const double *const b, double* x, uint32_t N_nod);

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
			 double *stress,
			 double *displacement, /* Output */
			 uint32_t N_force_steps,
			 bool *plastic_elements,
			 double *nodal_strain,
			 double *nodal_stress);

void get_stress_params(double *max_vm_stress, uint32_t *N_plastic_elem,
                       uint32_t *plastified_elem, double *stress, double yield_stress,
                       nb_plastified_analysis2D *elem_regime, uint32_t N_elem);

void find_new_plastic_elem(uint32_t N_elem,
                           nb_plastified_analysis2D *elem_regime,
                           double *stress, double yield_stress,
                           double *max_vm_stress,
                           uint32_t *plastified_elem,
                           uint32_t *N_plastic_elem);

char* strg_concat(char* print_name, char* problem_name, char* extension);

void interpolate_trg_strain_nodes_to_elemes(const nb_mesh2D_t *const part, double *strain, double *nodal_strain);
void interpolate_trg_stress_nodes_to_elemes(const nb_mesh2D_t *const part, double *stress, double *nodal_stress);
#endif // STATIC_PLASTICITY2D_H_INCLUDED
