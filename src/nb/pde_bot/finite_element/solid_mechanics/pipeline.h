#ifndef __NB_PDE_BOT_FINITE_ELEMENT_SOLID_MECHANICS_PIPELINE_H__
#define __NB_PDE_BOT_FINITE_ELEMENT_SOLID_MECHANICS_PIPELINE_H__

#include "nb/solver_bot.h"
#include "nb/geometric_bot.h"
#include "nb/pde_bot/material.h"
#include "nb/pde_bot/boundary_conditions/bcond.h"
#include "nb/pde_bot/finite_element/element.h"

bool pipeline_elem_is_enabled(const bool *elements_enabled, uint32_t id);
void pipeline_sum_gauss_point(const nb_fem_elem_t *elem, int gp_id,
			      double D[4], double density,
			      double thickness, double detJ,
			      double *dNi_dx, double *dNi_dy,
			      double fx, double fy,
			      double *Ke, double *Me, double *Fe);
void pipeline_add_to_global_system(const nb_fem_elem_t *elem, uint32_t id,
				   const nb_partition_t *part,
				   double *Ke, double *Me, double *Fe,
				   nb_sparse_t *K, double *M, double *F);
int pipeline_assemble_system
		(nb_sparse_t* K, double* M, double *F,
		 const nb_partition_t *const part,
		 const nb_fem_elem_t *const elemtype,
		 const nb_material_t *const material,
		 bool enable_self_weight,
		 double gravity[2],
		 nb_analysis2D_t analysis2D,
		 nb_analysis2D_params *params2D,
		 const bool* elements_enabled /* NULL to enable all */);

void pipeline_set_boundary_conditions(const nb_partition_t *part,
				      nb_sparse_t* K,
				      double* F, 
				      const nb_bcond_t *const bcond, 
				      double factor);
void pipeline_compute_strain(double *strain,
			     const nb_partition_t *const part,
			     double *displacement,
			     const nb_fem_elem_t *const elemtype);

void pipeline_compute_main_stress(double *stress, 
				  double *main_stress,
				  uint32_t N_elements,
				  const nb_fem_elem_t *const elemtype);

#endif
