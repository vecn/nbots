#ifndef __NB_PDE_BOT_FINITE_ELEMENT_SOLID_MECHANICS_PIPELINE_H__
#define __NB_PDE_BOT_FINITE_ELEMENT_SOLID_MECHANICS_PIPELINE_H__

#include "nb/eigen_bot.h"
#include "nb/geometric_bot.h"
#include "nb/pde_bot/material.h"
#include "nb/pde_bot/boundary_conditions/bcond.h"
#include "nb/pde_bot/finite_element/element.h"

bool pipeline_elem_is_enabled(const bool *elements_enabled, uint32_t id);
void pipeline_get_constitutive_matrix(double D[4], 
				      const vcn_fem_material_t *material,
				      nb_analysis2D_t analysis2D);
void pipeline_sum_gauss_point(const vcn_fem_elem_t *elem, int gp_id,
			      double D[4], double density,
			      double thickness, double detJ,
			      double *dNi_dx, double *dNi_dy,
			      double fx, double fy,
			      double *Ke, double *Me, double *Fe);
void pipeline_add_to_global_system(const vcn_fem_elem_t *elem, uint32_t id,
				   const vcn_msh3trg_t *mesh,
				   double *Ke, double *Me, double *Fe,
				   vcn_sparse_t *K, double *M, double *F);
int pipeline_assemble_system
		(vcn_sparse_t* K, double* M, double *F,
		 const vcn_msh3trg_t *const mesh,
		 const vcn_fem_elem_t *const elemtype,
		 const vcn_fem_material_t *const material,
		 bool enable_self_weight,
		 double gravity[2],
		 nb_analysis2D_t analysis2D,
		 nb_analysis2D_params *params2D,
		 const bool* elements_enabled /* NULL to enable all */);

void pipeline_set_boundary_conditions(const vcn_msh3trg_t *msh3trg,
				      vcn_sparse_t* K,
				      double* F, 
				      const nb_bcond_t *const bcond, 
				      double factor);
void pipeline_compute_strain(double *strain,
			     const vcn_msh3trg_t *const mesh,
			     double *displacement,
			     const vcn_fem_elem_t *const elemtype,
			     nb_analysis2D_t analysis2D,
			     const vcn_fem_material_t *const material);

void pipeline_compute_main_stress(double *stress, 
				  double *main_stress,
				  uint32_t N_elements,
				  const vcn_fem_elem_t *const elemtype);

#endif
