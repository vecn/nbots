#ifndef __VCN_PDE_BOT_FINITE_ELEMENT_SOLID_MECHANICS_PIPELINE_H__
#define __VCN_PDE_BOT_FINITE_ELEMENT_SOLID_MECHANICS_PIPELINE_H__

void pipeline_assemble_system
		(vcn_sparse_t* K, double* M, double *F,
		 const vcn_msh3trg_t *const mesh,
		 const vcn_fem_elem_t *const elemtype,
		 const vcn_fem_material_t *const material,
		 bool enable_self_weight,
		 double gravity[2],
		 bool enable_plane_stress,
		 double thickness,
		 bool* elements_enabled /* NULL to enable all */);

void pipeline_set_boundary_conditions(const vcn_msh3trg_t *msh3trg,
				      vcn_sparse_t* K,
				      double* F, 
				      const vcn_bcond_t *const bmeshcond, 
				      double thickness,
				      double factor);
void pipeline_compute_strain
			(double *strain,
			 const vcn_msh3trg_t *const mesh,
			 double *displacement,
			 const vcn_fem_elem_t *const elemtype,
			 bool enable_plane_stress,
			 const vcn_fem_material_t *const material);

void pipeline_compute_main_stress(double *stress, 
				  double *main_stress,
				  uint32_t N_elements,
				  const vcn_fem_elem_t *const elemtype);

void pipeline_compute_error_on_elements
			(double* error,
			 const vcn_msh3trg_t *const mesh,
			 double *displacement,
			 double* strain,
			 const vcn_fem_elem_t *const elemtype);

#endif
