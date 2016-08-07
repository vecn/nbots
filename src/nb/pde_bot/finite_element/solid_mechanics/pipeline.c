#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <alloca.h>

#include "nb/math_bot.h"
#include "nb/cfreader_cat.h"
#include "nb/eigen_bot.h"
#include "nb/container_bot.h"
#include "nb/graph_bot.h"
#include "nb/pde_bot/material.h"
#include "nb/pde_bot/boundary_conditions/bcond.h"
#include "nb/pde_bot/boundary_conditions/bcond_iter.h"
#include "nb/pde_bot/finite_element/element.h"
#include "nb/pde_bot/finite_element/gaussp_to_nodes.h"
#include "nb/pde_bot/finite_element/solid_mechanics/static_elasticity2D.h"
  
#include "../element_struct.h"
#include "pipeline.h"

#define POW2(a) ((a)*(a))

static bool elem_is_enabled(const bool *elements_enabled, uint32_t id);
static int assemble_element(const vcn_fem_elem_t *elem, uint32_t id,
			    const vcn_msh3trg_t *mesh,
			    const vcn_fem_material_t *material,
			    bool is_enabled,
			    nb_analysis2D_t analysis2D,
			    nb_analysis2D_params *params2D,
			    bool enable_self_weight,
			    double gravity[2],
			    vcn_sparse_t *K, double *M, double *F);
static void get_constitutive_matrix(double D[4], 
				    const vcn_fem_material_t *material,
				    nb_analysis2D_t analysis2D);
static void set_plane_stress(double D[4], double E, double v);
static void set_plane_strain(double D[4], double E, double v);
static int integrate_elemental_system
		       	(const vcn_fem_elem_t *elem, uint32_t id,
			 double D[4], double density, double gravity[2],
			 const vcn_msh3trg_t *mesh,
			 nb_analysis2D_params *params2D,
			 bool enable_self_weight,
			 double *Ke, double *Me, double *Fe);
static void set_bcond_neumann_sgm(const vcn_msh3trg_t *msh3trg,
				  vcn_sparse_t* K, double* F, 
				  const nb_bcond_t *const bcond, 
				  double factor);
static void set_bcond_neumann_vtx(const vcn_msh3trg_t *msh3trg,
				  vcn_sparse_t* K, double* F, 
				  const nb_bcond_t *const bcond, 
				  double factor);
static void set_bcond_dirichlet_sgm(const vcn_msh3trg_t *msh3trg,
				    vcn_sparse_t* K, double* F, 
				    const nb_bcond_t *const bcond, 
				    double factor);
static void set_bcond_dirichlet_vtx(const vcn_msh3trg_t *msh3trg,
				    vcn_sparse_t* K, double* F, 
				    const nb_bcond_t *const bcond, 
				    double factor);
int pipeline_assemble_system
		(vcn_sparse_t* K, double* M, double *F,
		 const vcn_msh3trg_t *const mesh,
		 const vcn_fem_elem_t *const elem,
		 const vcn_fem_material_t *const material,
		 bool enable_self_weight,
		 double gravity[2],
		 nb_analysis2D_t analysis2D,
		 nb_analysis2D_params *params2D,
		 const bool* elements_enabled /* NULL to enable all */)
{
	int status = 1;
	uint32_t N_elements = mesh->N_triangles;
	vcn_sparse_reset(K);
	if (NULL != M)
		memset(M, 0, vcn_sparse_get_size(K) * sizeof(*M));
	memset(F, 0, vcn_sparse_get_size(K) * sizeof(*F));

	for (uint32_t k = 0; k < N_elements; k++) {
		bool is_enabled = elem_is_enabled(elements_enabled, k);
		int status_element =
			assemble_element(elem, k, mesh, material, is_enabled,
					 analysis2D, params2D,
					 enable_self_weight, gravity,
					 K, M, F);
		if (0 != status_element)
			goto EXIT;
	}
	status = 0;
EXIT:
	return status;
}

static bool elem_is_enabled(const bool *elements_enabled, uint32_t id)
{
	bool is_enabled = true;
	if (NULL != elements_enabled)
		is_enabled = elements_enabled[id];
	return is_enabled;
}

static int assemble_element(const vcn_fem_elem_t *elem, uint32_t id,
			    const vcn_msh3trg_t *mesh,
			    const vcn_fem_material_t *material,
			    bool is_enabled,
			    nb_analysis2D_t analysis2D,
			    nb_analysis2D_params *params2D,
			    bool enable_self_weight,
			    double gravity[2],
			    vcn_sparse_t *K, double *M, double *F)
{
	int status = 1;
    
	double D[4] = {1e-6, 1e-6, 1e-6, 1e-6};
	double density = 1e-6;
	if (is_enabled) {
		get_constitutive_matrix(D, material, analysis2D);
		density = vcn_fem_material_get_density(material);
	}

	/* Allocate elemental Stiffness Matrix and Force Vector */
	double* Ke = calloc(4 * POW2(elem->N_nodes), sizeof(*Ke));
	double* Me = NULL;
	if(M != NULL)
		Me = calloc(2 * elem->N_nodes, sizeof(*Me));
	double* Fe = calloc(2 * elem->N_nodes, sizeof(*Fe));
	
	int status_elem = integrate_elemental_system(elem, id, D, density,
						     gravity, mesh, params2D,
						     enable_self_weight,
						     Ke, Me, Fe);
	if (0 != status_elem)
		goto CLEANUP;
	
	/* Add to global stiffness matrix */
	for (uint32_t i1 = 0; i1 < elem->N_nodes; i1++) {
		for (uint32_t i2 = 0; i2 < elem->N_nodes; i2++) {
			for (uint8_t j1=0; j1 < 2; j1++) {
				for (uint8_t j2=0; j2 < 2; j2++) {
					vcn_sparse_add
						(K, mesh->vertices_forming_triangles[id*3+i1]*2 + j1,
						 mesh->vertices_forming_triangles[id*3+i2]*2 + j2,
						 Ke[(i1*2+j1)*(2*elem->N_nodes) + (i2*2+j2)]);
				}
			}
		}
		/* Add to global mass matrix */
		if (NULL != M) {
			M[mesh->vertices_forming_triangles[id*3+i1] * 2] += Me[i1 * 2];
			M[mesh->vertices_forming_triangles[id*3+i1]*2+1] += Me[i1*2+1];
		}
		/* Add to global forces vector */
		F[mesh->vertices_forming_triangles[id*3+i1] * 2] += Fe[i1 * 2];
		F[mesh->vertices_forming_triangles[id*3+i1]*2+1] += Fe[i1*2+1];
	}
	status = 0;
CLEANUP:
	free(Ke);
	if (NULL != M) 
		free(Me);
	free(Fe);
	return status;
}

static void get_constitutive_matrix(double D[4], 
				    const vcn_fem_material_t *material,
				    nb_analysis2D_t analysis2D)
{
	double E = vcn_fem_material_get_elasticity_module(material);
	double v = vcn_fem_material_get_poisson_module(material);
	switch (analysis2D) {
	case NB_PLANE_STRESS:
		set_plane_stress(D, E, v);
	case NB_PLANE_STRAIN:
		set_plane_strain(D, E, v);
	default:
		set_plane_stress(D, E, v);
	}
}

static void set_plane_stress(double D[4], double E, double v)
{
	D[0] = E / (1.0 - POW2(v));
	D[1] = v * D[0];
	D[2] = D[0];
	D[3] = E / (2.0 * (1.0 + v));
}

static void set_plane_strain(double D[4], double E, double v)
{
	D[0] = (E * (1.0 - v)) / ((1.0 + v) * (1.0 - 2 * v));
	D[1] = (v * D[0]) / (1.0 - v);
	D[2] = D[0];
	D[3] = E / (2.0 * (1.0 + v));
}

static int integrate_elemental_system
		       	(const vcn_fem_elem_t *elem, uint32_t id,
			 double D[4], double density, double gravity[2],
			 const vcn_msh3trg_t *mesh,
			 nb_analysis2D_params *params2D,
			 bool enable_self_weight,
			 double *Ke, double *Me, double *Fe)
{
	int status = 1;

	double fx = 0.0;
	double fy = 0.0;
	if (enable_self_weight) {
		fx = gravity[0] * density;
		fy = gravity[1] * density;
	}
    
	/* Allocate Cartesian derivatives for each Gauss Point */
	double *dNi_dx = malloc(elem->N_nodes * sizeof(*dNi_dx));
	double *dNi_dy = malloc(elem->N_nodes * sizeof(*dNi_dy));
	for (uint32_t j = 0; j < elem->N_Gauss_points; j++) {      
		/* Compute Jacobian derivatives */
		double dx_dpsi = 0.0;
		double dy_dpsi = 0.0;
		double dx_deta = 0.0;
		double dy_deta = 0.0;
		for (uint32_t i = 0; i < elem->N_nodes; i++) {
			uint32_t inode = mesh->vertices_forming_triangles[id * elem->N_nodes + i];
			double xi = mesh->vertices[inode * 2];
			double yi = mesh->vertices[inode*2+1];
			dx_dpsi += 
				elem->dNi_dpsi[i](elem->psi[j], elem->eta[j]) * xi;
			dx_deta +=
				elem->dNi_deta[i](elem->psi[j], elem->eta[j]) * xi;
			dy_dpsi += 
				elem->dNi_dpsi[i](elem->psi[j], elem->eta[j]) * yi;
			dy_deta +=
				elem->dNi_deta[i](elem->psi[j], elem->eta[j]) * yi;
		}

		/* Compute Jacobian inverse and determinant */
		double detJ = dx_dpsi * dy_deta - dy_dpsi * dx_deta;
		/* Check if the element is distorted */
		if(detJ < 0)
			goto EXIT;
      
		double Jinv[4];
		Jinv[0] =  dy_deta/detJ;
		Jinv[1] = -dy_dpsi/detJ;
		Jinv[2] = -dx_deta/detJ;
		Jinv[3] =  dx_dpsi/detJ;
      
		/* Compute Shape functions derivatives in cartesian space */ 
		for (uint32_t i = 0; i < elem->N_nodes; i++) {
			dNi_dx[i] = 
				Jinv[0]*elem->dNi_dpsi[i](elem->psi[j], elem->eta[j]) + 
				Jinv[1]*elem->dNi_deta[i](elem->psi[j], elem->eta[j]);
			dNi_dy[i] = 
				Jinv[2]*elem->dNi_dpsi[i](elem->psi[j], elem->eta[j]) + 
				Jinv[3]*elem->dNi_deta[i](elem->psi[j], elem->eta[j]);
		}
		/* Compute elemental stiffness matrix
		 *       _            _         _           _        _  _
		 *      |dNi/dx        |       | d11 d12     |      |  |
		 * Bi = |       dNi/dy |   D = | d21 d22     |  K = |  | B'DB dx dy
		 *      |dNi/dy dNi/dx_|       |_        d33_|     _| _|

		 * B  = [B1 B2...Bi... Belem->N_nodes]
		 *
		 */
      
		for (uint32_t i1 = 0; i1 < elem->N_nodes; i1++) {
			for (uint32_t i2 = 0; i2 < elem->N_nodes; i2++) {
				/*  Integrating elemental siffness matrix */
				Ke[(i1 * 2)*(2*elem->N_nodes) + (i2 * 2)] += 
					(dNi_dx[i1]*dNi_dx[i2]*D[0] + dNi_dy[i1]*dNi_dy[i2]*D[3]) *
					detJ * params2D->thickness * elem->gp_weight[j];
				Ke[(i1 * 2)*(2*elem->N_nodes) + (i2*2+1)] +=
					(dNi_dx[i1]*dNi_dy[i2]*D[1] + dNi_dy[i1]*dNi_dx[i2]*D[3]) *
					detJ * params2D->thickness * elem->gp_weight[j];
				Ke[(i1*2+1)*(2*elem->N_nodes) + (i2 * 2)] +=
					(dNi_dy[i1]*dNi_dx[i2]*D[1] + dNi_dx[i1]*dNi_dy[i2]*D[3]) *
					detJ * params2D->thickness * elem->gp_weight[j];
				Ke[(i1*2+1)*(2*elem->N_nodes) + (i2*2+1)] +=
					(dNi_dy[i1]*dNi_dy[i2]*D[2] + dNi_dx[i1]*dNi_dx[i2]*D[3]) *
					detJ * params2D->thickness * elem->gp_weight[j];
			}
			/* Calculate shape function of the i-th node at the j-th gauss point */
			double Ni_eval = 
				elem->Ni[i1](elem->psi[j], elem->eta[j]);
			/*  Integrating elemental mass matrix */
			if (NULL != Me) {
				/* OPPORTUNITY: Allocate just one for each component */
				Me[i1 * 2] += Ni_eval * Ni_eval * density *
					detJ * params2D->thickness * elem->gp_weight[j];
				Me[i1*2+1] += Ni_eval * Ni_eval * density * 
					detJ * params2D->thickness * elem->gp_weight[j];
			}
			/* Compute elemental forces vector */
			/* OPPORTUNITY: Allocate just one for each component */
			Fe[i1 * 2] += Ni_eval * fx * detJ * 
				params2D->thickness * elem->gp_weight[j];
			Fe[i1*2+1] += Ni_eval * fy * detJ * 
				params2D->thickness * elem->gp_weight[j];
		}
	}
	status = 0;
EXIT:
	free(dNi_dx);
	free(dNi_dy);
	return status;
}

void pipeline_set_boundary_conditions(const vcn_msh3trg_t *msh3trg,
				      vcn_sparse_t* K, double* F, 
				      const nb_bcond_t *const bcond,
				      double factor)
{
	set_bcond_neumann_sgm(msh3trg, K, F, bcond, factor);
	set_bcond_neumann_vtx(msh3trg, K, F, bcond, factor);
	set_bcond_dirichlet_sgm(msh3trg, K, F, bcond, factor);
	set_bcond_dirichlet_vtx(msh3trg, K, F, bcond, factor);
}

static void set_bcond_neumann_sgm(const vcn_msh3trg_t *msh3trg,
				  vcn_sparse_t* K,
				  double* F, 
				  const nb_bcond_t *const bcond, 
				  double factor)
{
	uint8_t N_dof = nb_bcond_get_N_dof(bcond);
	uint16_t size = nb_bcond_iter_get_memsize();
	nb_bcond_iter_t *iter = alloca(size);
	nb_bcond_iter_init(iter);
	nb_bcond_iter_set_conditions(iter, bcond, NB_NEUMANN,
				     NB_BC_ON_SEGMENT);
	while (nb_bcond_iter_has_more(iter)) {
		nb_bcond_iter_go_next(iter);
		uint32_t model_id = nb_bcond_iter_get_id(iter);
		uint32_t N = msh3trg->N_vtx_x_inputsgm[model_id];
		for (uint32_t i = 0; i < N; i++) {
			uint32_t mesh_id = 
				msh3trg->meshvtx_x_inputsgm[model_id][i];
			for (uint32_t j = 0; j < N_dof; j++) {
				bool mask = nb_bcond_iter_get_mask(iter, j);
				double val = nb_bcond_iter_get_val(iter, j);
				if (mask) {
					uint32_t mtx_id = mesh_id * N_dof + j;
					F[mtx_id] += factor * (val / N);
				}
			}
		}
	}
	nb_bcond_iter_finish(iter);
}

static void set_bcond_neumann_vtx(const vcn_msh3trg_t *msh3trg,
				  vcn_sparse_t* K,
				  double* F, 
				  const nb_bcond_t *const bcond,
				  double factor)
{
	uint8_t N_dof = nb_bcond_get_N_dof(bcond);
	uint16_t size = nb_bcond_iter_get_memsize();
	nb_bcond_iter_t *iter = alloca(size);
	nb_bcond_iter_init(iter);
	nb_bcond_iter_set_conditions(iter, bcond, NB_NEUMANN,
				     NB_BC_ON_POINT);
	while (nb_bcond_iter_has_more(iter)) {
		nb_bcond_iter_go_next(iter);
		uint32_t model_id = nb_bcond_iter_get_id(iter);
		uint32_t mesh_id = msh3trg->input_vertices[model_id];
		for (uint32_t j = 0; j < N_dof; j++) {
			bool mask = nb_bcond_iter_get_mask(iter, j);
			double val = nb_bcond_iter_get_val(iter, j);
			if (mask) {
				uint32_t mtx_id = mesh_id * N_dof + j;
				F[mtx_id] += factor * val;
			}
		}
	}
	nb_bcond_iter_finish(iter);
}

static void set_bcond_dirichlet_sgm(const vcn_msh3trg_t *msh3trg,
				    vcn_sparse_t* K,
				    double* F, 
				    const nb_bcond_t *const bcond,
				    double factor)
{
	uint8_t N_dof = nb_bcond_get_N_dof(bcond);
	uint16_t size = nb_bcond_iter_get_memsize();
	nb_bcond_iter_t *iter = alloca(size);
	nb_bcond_iter_init(iter);
	nb_bcond_iter_set_conditions(iter, bcond, NB_DIRICHLET,
				     NB_BC_ON_SEGMENT);
	while (nb_bcond_iter_has_more(iter)) {
		nb_bcond_iter_go_next(iter);
		uint32_t model_id = nb_bcond_iter_get_id(iter);
		uint32_t N = msh3trg->N_vtx_x_inputsgm[model_id];
		for (uint32_t i = 0; i < N; i++) {
			uint32_t mesh_id = 
				msh3trg->meshvtx_x_inputsgm[model_id][i];
			for (uint32_t j = 0; j < N_dof; j++) {
				bool mask = nb_bcond_iter_get_mask(iter, j);
				double val = nb_bcond_iter_get_val(iter, j);
				if (mask) {
					uint32_t mtx_id = mesh_id * N_dof + j;
					vcn_sparse_set_Dirichlet_condition
						(K, F, mtx_id, factor * val);
				}
			}
		}
	}
	nb_bcond_iter_finish(iter);
}

static void set_bcond_dirichlet_vtx(const vcn_msh3trg_t *msh3trg,
				    vcn_sparse_t* K,
				    double* F, 
				    const nb_bcond_t *const bcond,
				    double factor)
{
	uint8_t N_dof = nb_bcond_get_N_dof(bcond);
	uint16_t size = nb_bcond_iter_get_memsize();
	nb_bcond_iter_t *iter = alloca(size);
	nb_bcond_iter_init(iter);
	nb_bcond_iter_set_conditions(iter, bcond, NB_DIRICHLET,
				     NB_BC_ON_POINT);
	while (nb_bcond_iter_has_more(iter)) {
		nb_bcond_iter_go_next(iter);
		uint32_t model_id = nb_bcond_iter_get_id(iter);
		uint32_t mesh_id = msh3trg->input_vertices[model_id];
		for (uint32_t j = 0; j < N_dof; j++) {
			bool mask = nb_bcond_iter_get_mask(iter, j);
			double val = nb_bcond_iter_get_val(iter, j);
			if (mask) {
				uint32_t mtx_id = mesh_id * N_dof + j;
				vcn_sparse_set_Dirichlet_condition
					(K, F, mtx_id, factor * val);
			}
		}
	}
	nb_bcond_iter_finish(iter);
}

void pipeline_compute_strain(double *strain,
			     const vcn_msh3trg_t *const mesh,
			     double *displacement,
			     const vcn_fem_elem_t *const elem,
			     nb_analysis2D_t analysis2D,
			     const vcn_fem_material_t *const material)
{
	double *vertices = (double*) mesh->vertices;
	uint32_t N_elements = mesh->N_triangles;
	uint32_t *connectivity_mtx = (uint32_t*) mesh->vertices_forming_triangles;

	/* Initialize strains */
	memset(strain, 0, 3 * elem->N_Gauss_points * N_elements * sizeof(*strain));

	/* Iterate over elements to compute strain and stress at nodes */
	for (uint32_t k = 0; k < N_elements; k++) {
		double* dNi_dx = (double*)malloc(elem->N_nodes*sizeof(double));
		double* dNi_dy = (double*)malloc(elem->N_nodes*sizeof(double));

		/* Integrate domain */
		for(uint32_t j=0; j < elem->N_Gauss_points; j++){
			uint32_t idx = k * elem->N_Gauss_points + j;

			/* Compute Jacobian derivatives */
			double dx_dpsi = 0.0;
			double dy_dpsi = 0.0;
			double dx_deta = 0.0;
			double dy_deta = 0.0;

			for(uint32_t i=0; i < elem->N_nodes; i++){
				uint32_t inode = connectivity_mtx[k * elem->N_nodes + i];
				double xi = vertices[inode * 2];
				double yi = vertices[inode*2+1];
				dx_dpsi += 
					elem->dNi_dpsi[i](elem->psi[j], elem->eta[j]) * xi;
				dx_deta += 
					elem->dNi_deta[i](elem->psi[j], elem->eta[j]) * xi;
				dy_dpsi +=
					elem->dNi_dpsi[i](elem->psi[j], elem->eta[j]) * yi;
				dy_deta +=
					elem->dNi_deta[i](elem->psi[j], elem->eta[j]) * yi;
			}
      
			/* Compute Jacobian inverse and determinant */
			double detJ = dx_dpsi*dy_deta - dy_dpsi*dx_deta;
			double Jinv[4];
			Jinv[0] =  dy_deta/detJ;
			Jinv[1] = -dy_dpsi/detJ;
			Jinv[2] = -dx_deta/detJ;
			Jinv[3] =  dx_dpsi/detJ;
      
			/* Compute Shape functions derivatives in cartesian space */ 
			for(uint32_t i=0; i < elem->N_nodes; i++){
				dNi_dx[i] = 
					Jinv[0]*elem->dNi_dpsi[i](elem->psi[j], elem->eta[j]) + 
					Jinv[1]*elem->dNi_deta[i](elem->psi[j], elem->eta[j]);
				dNi_dy[i] = 
					Jinv[2]*elem->dNi_dpsi[i](elem->psi[j], elem->eta[j]) + 
					Jinv[3]*elem->dNi_deta[i](elem->psi[j], elem->eta[j]);
			}

			/* Compute Strain at Gauss Point */
			for (uint32_t i=0; i < elem->N_nodes; i++) {
				uint32_t inode = connectivity_mtx[k * elem->N_nodes + i];
				strain[idx * 3] += dNi_dx[i] * displacement[inode * 2];
				strain[idx*3+1] += dNi_dy[i] * displacement[inode*2+1];
				strain[idx*3+2] += (dNi_dy[i] * displacement[inode * 2] +
						    dNi_dx[i] * displacement[inode*2+1]);
			}
      		}
		free(dNi_dx);
		free(dNi_dy);
	}
}

void pipeline_compute_main_stress(double *stress, 
				  double *main_stress,
				  uint32_t N_elements,
				  const vcn_fem_elem_t *const elem)
{
	for(uint32_t i=0; i < N_elements; i++){
		for(uint32_t j=0; j < elem->N_Gauss_points; j++){
			uint32_t idx = i*elem->N_Gauss_points + j;
			double sigma_avg = 0.5 * (stress[idx * 3] + stress[idx*3+1]);
			double R = sqrt(0.25 * 
					POW2(stress[idx * 3] - stress[idx*3+1]) +
					POW2(stress[idx*3+2]));
			main_stress[idx * 2] = sigma_avg + R;
			main_stress[idx*2+1] = sigma_avg - R;
		}
	}
}

void pipeline_compute_error_on_elements
			(double* error,
			 const vcn_msh3trg_t *const mesh,
			 double *displacement,
			 double* strain,
			 const vcn_fem_elem_t *const elem)
/* Compute elemental error based on strain field */
{
	uint32_t N_vertices = mesh->N_vertices;
	uint32_t N_elements = mesh->N_triangles;


	uint32_t N_gp = elem->N_Gauss_points;
	memset(error, 0, N_elements * N_gp * sizeof(double));

	/* Interpolate strain to nodes */
	double* strain_on_nodes = (double*)calloc(3 * N_vertices, sizeof(double));
	vcn_fem_interpolate_from_Gauss_points_to_nodes(mesh, elem,
						       3, strain, strain_on_nodes);

	uint32_t *connectivity_mtx = (uint32_t*) mesh->vertices_forming_triangles;

	/* Return values to Gauss Points to compute the error in the element */
	for(uint32_t i=0; i < N_elements; i++){
		error[i] = 0;
		for(uint32_t j=0; j < N_gp; j++){
			double strain_gp[3];
			memset(strain_gp, 0, 3*sizeof(double));
			for(uint32_t k=0; k < elem->N_nodes; k++){
				strain_gp[0] += elem->Ni[k](elem->psi[j], elem->eta[j]) * 
					strain_on_nodes[connectivity_mtx[i*3+k] * 3];
				strain_gp[1] += elem->Ni[k](elem->psi[j], elem->eta[j]) * 
					strain_on_nodes[connectivity_mtx[i*3+k]*3+1];
				strain_gp[2] += elem->Ni[k](elem->psi[j], elem->eta[j]) * 
					strain_on_nodes[connectivity_mtx[i*3+k]*3+2];
			}
			uint32_t idx = i * N_gp + j;
			error[idx] += 
				sqrt(POW2(strain[idx * 3] - strain_gp[0]) +
				     POW2(strain[idx*3+1] - strain_gp[1]) +
				     POW2(strain[idx*3+2] - strain_gp[2]));
		}
	}

	/* Free memory */
	free(strain_on_nodes);
}
