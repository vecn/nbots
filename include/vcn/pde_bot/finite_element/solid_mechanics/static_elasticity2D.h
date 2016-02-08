#ifndef __VCN_PDE_BOT_FINITE_ELEMENT_SOLID_MECHANICS_STATIC_ELASTICITY2D_H__
#define __VCN_PDE_BOT_FINITE_ELEMENT_SOLID_MECHANICS_STATIC_ELASTICITY2D_H__

#include <stdint.h>
#include <stdbool.h>
#include "vcn/eigen_bot.h"
#include "vcn/geometric_bot.h"

#include "vcn/pde_bot/material.h"
#include "vcn/pde_bot/boundary_conditions.h"
#include "vcn/pde_bot/finite_element/element.h"

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
			 double* displacement, /* Output */
			 double* strain,       /* Output */
			 const char* logfile /* NULL if not required */);

void vcn_fem_compute_stress_from_strain
			(uint32_t N_elements,
			 uint32_t* elements_connectivity_matrix, 
			 const vcn_fem_elem_t *const elemtype,
			 const vcn_fem_material_t *const material,
			 bool enable_plane_stress,
			 double* strain,
			 uint32_t omp_parallel_threads,
			 bool* elements_enabled /* NULL to enable all */,
			 double* stress /* Output */);

#endif
