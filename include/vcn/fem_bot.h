/******************************************************************************
 *   FEM Bot: Finit element Method Bot                                        *
 *   2011-2015 Victor Eduardo Cardoso Nungaray                                *
 *   Twitter: @victore_cardoso                                                *
 *   email: victorc@cimat.mx                                                  *
 ******************************************************************************/

/**
 * @file fem_bot.h
 * @brief FEM bot will assemble the linear systems of well stablished FE 
 * problems, such as solid mechanics and heat diffusion.
 * @author Victor Eduardo Cardoso Nungaray
 * @n victorc@@cimat.mx
 * @n @@victore_cardoso
 * @date 10 August 2015
 */

#ifndef __VCN_FEM_BOT_H__
#define __VCN_FEM_BOT_H__

#include <stdint.h>
#include <stdbool.h>
#include "vcn/eigen_bot.h"
#include "vcn/geometric_bot.h"

#ifdef __cplusplus
extern "C" {
#endif
  
  /************ Types definition ****************/ 
  typedef struct vcn_fem_elem_s vcn_fem_elem_t;
  typedef struct vcn_fem_material_s vcn_fem_material_t;
  typedef struct vcn_fem_implicit_s vcn_fem_implicit_t;

  typedef struct {
    /* Boundary conditions */
    uint8_t N_dof; /* Degrees of freedom */
    /* Dirichlet conditions on vertices */
    uint32_t N_Dirichlet_on_vtx;
    uint32_t *Dirichlet_on_vtx_idx;
    bool *Dirichlet_on_vtx_dof_mask;
    double *Dirichlet_on_vtx_val;
    /* Neuman conditions on vertices */
    uint32_t N_Neuman_on_vtx;
    uint32_t *Neuman_on_vtx_idx;
    bool *Neuman_on_vtx_dof_mask;
    double *Neuman_on_vtx_val;
    /* Dirichlet conditions on segments */
    uint32_t N_Dirichlet_on_sgm;
    uint32_t *Dirichlet_on_sgm_idx;
    bool*Dirichlet_on_sgm_dof_mask;
    double *Dirichlet_on_sgm_val;
    /* Neuman conditions on segments */
    uint32_t N_Neuman_on_sgm;
    uint32_t *Neuman_on_sgm_idx;
    bool *Neuman_on_sgm_dof_mask;
    double *Neuman_on_sgm_val;
  }vcn_bcond_t;

  typedef struct{
    /* General information */
    char filename[255];
    char author[255];
    char project_name[255];
    bool enable_double_precision;
  }vcn_fem_output_t;

  /************ Functions definition ************/
  vcn_fem_implicit_t* vcn_fem_implicit_create();
  void vcn_fem_implicit_destroy(vcn_fem_implicit_t* isparams);

  void vcn_fem_implicit_set_N_steps(vcn_fem_implicit_t* isparams,
					       uint32_t N_time_steps);
  void vcn_fem_implicit_set_N_max_iter(vcn_fem_implicit_t* isparams,
					     uint32_t N_max_iter);
  void vcn_fem_implicit_set_N_max_iter_without_enhance
  (vcn_fem_implicit_t* isparams, uint32_t N_max_iter);
  void vcn_fem_implicit_set_residual_tolerance
       (vcn_fem_implicit_t* isparams, double penergy_tol);

  uint32_t vcn_fem_implicit_get_N_steps(vcn_fem_implicit_t* isparams);
  uint32_t vcn_fem_implicit_get_N_max_iter(vcn_fem_implicit_t* isparams);
  uint32_t vcn_fem_implicit_get_N_max_iter_without_enhance
  (vcn_fem_implicit_t* isparams);
  double vcn_fem_implicit_get_residual_tolerance
         (vcn_fem_implicit_t* isparams);
  
  vcn_fem_material_t* vcn_fem_material_create();
  void vcn_fem_material_destroy(vcn_fem_material_t* mat);

  void vcn_fem_material_set_poisson_module(vcn_fem_material_t* mat, double pmodule);
  void vcn_fem_material_set_elasticity_module(vcn_fem_material_t* mat, double emodule);
  void vcn_fem_material_set_density(vcn_fem_material_t* mat, double density);
  void vcn_fem_material_set_fracture_energy(vcn_fem_material_t* mat, double frac_energy);
  void vcn_fem_material_set_traction_limit_stress(vcn_fem_material_t* mat, double max_stress);
  void vcn_fem_material_set_compression_limit_stress(vcn_fem_material_t* mat, double max_stress);

  double vcn_fem_material_get_poisson_module(const vcn_fem_material_t *const mat);
  double vcn_fem_material_get_elasticity_module(const vcn_fem_material_t *const mat);
  double vcn_fem_material_get_density(const vcn_fem_material_t *const mat);
  double vcn_fem_material_get_fracture_energy(const vcn_fem_material_t *const mat);
  double vcn_fem_material_get_traction_limit_stress(const vcn_fem_material_t *const mat);
  double vcn_fem_material_get_compression_limit_stress(const vcn_fem_material_t *const mat);
  void* vcn_fem_material_get_damage_function(const vcn_fem_material_t *const mat);

  int vcn_fem_material_verify(vcn_fem_material_t* material);

  vcn_bcond_t* vcn_fem_bcond_create();
  vcn_bcond_t* vcn_fem_bcond_clone(const vcn_bcond_t *const bconditions);
  vcn_bcond_t* vcn_fem_bcond_read(const char* filename);
  void vcn_fem_bcond_save(const vcn_bcond_t *const bconditions, 
			 const char* filename);
  void vcn_fem_bcond_copy(vcn_bcond_t* cpy,
			const vcn_bcond_t *const src);
  void vcn_fem_bcond_clear(vcn_bcond_t* bconditions);
  void vcn_fem_bcond_destroy(vcn_bcond_t* bconditions);
  void vcn_fem_bcond_printf(const vcn_bcond_t* const bconditions);

  vcn_bcond_t* vcn_fem_bcond_create_from_model_to_mesh
                        (const vcn_msh3trg_t *const msh3trg,
			 const vcn_bcond_t *const bconditions);

  vcn_fem_elem_t* vcn_fem_elem_create_triangle();
  void vcn_fem_elem_destroy(vcn_fem_elem_t* elemtype);
  uint32_t vcn_fem_elem_get_N_nodes(const vcn_fem_elem_t *const elemtype);
  void* vcn_fem_elem_get_ith_shape_function
                        (const vcn_fem_elem_t *const elemtype, uint32_t i);
  uint32_t vcn_fem_elem_get_closest_Gauss_Point_to_the_ith_node
               (const vcn_fem_elem_t *const elemtype, uint32_t i);

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
   const char* logfile /* NULL: Does not create an output logifle */);

  void vcn_fem_compute_2D_Non_Linear_Solid_Mechanics
  (const vcn_msh3trg_t *const mesh,
   const vcn_fem_elem_t *const elemtype,
   const vcn_fem_material_t *const material,
   const vcn_bcond_t *const bmeshcond,
   bool enable_self_weight,
   double gravity[2],
   bool enable_Cholesky_solver,
   bool enable_plane_stress,
   double thickness,
   vcn_fem_implicit_t* params,
   bool restore_computation, /* Restore computation after crash */
   vcn_fem_output_t* output,
   const char* logfile);

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

  void vcn_fem_interpolate_from_Gauss_points_to_nodes
  (const vcn_msh3trg_t *const mesh,
   const vcn_fem_elem_t *const elemtype,
   uint32_t N_components,
   double* values_on_GP_from_elements,
   double* values_interpolated_on_nodes /* Output */);

#ifdef __cplusplus
}
#endif

#endif
