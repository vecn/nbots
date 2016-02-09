#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "vcn/math_bot.h"
#include "vcn/cfreader_cat.h"
#include "vcn/eigen_bot.h"
#include "vcn/container_bot.h"
#include "vcn/graph_bot.h"
#include "vcn/pde_bot/material.h"
#include "vcn/pde_bot/boundary_conditions.h"
#include "vcn/pde_bot/finite_element/element.h"
#include "vcn/pde_bot/finite_element/gaussp_to_nodes.h"
#include "vcn/pde_bot/finite_element/solid_mechanics/static_elasticity2D.h"
  
#include "../element_struct.h"
#include "pipeline.h"

#define POW2(a) ((a)*(a))

void pipeline_assemble_system
		(vcn_sparse_t* K, double* M, double *F,
		 const vcn_msh3trg_t *const mesh,
		 const vcn_fem_elem_t *const elemtype,
		 const vcn_fem_material_t *const material,
		 bool enable_self_weight,
		 double gravity[2],
		 bool enable_plane_stress,
		 double thickness,
		 bool* elements_enabled /* NULL to enable all */)
{
	uint32_t N_elements = mesh->N_triangles;

	vcn_sparse_reset(K);
	if (NULL != M)
		memset(M, 0, vcn_sparse_get_size(K) * sizeof(double));
	memset(F, 0, vcn_sparse_get_size(K) * sizeof(double));
	/* Allocate elemental Stiffness Matrix and Force Vector */
	double* Ke = (double*)malloc(4 * elemtype->N_nodes * elemtype->N_nodes *
				     sizeof(double));
	double* Me = NULL;
	if(M != NULL)
		Me = (double*)malloc(2 * elemtype->N_nodes*sizeof(double));
	double* Fe = (double*)
		malloc(2 * elemtype->N_nodes*sizeof(double));

	if (!enable_plane_stress)
		thickness = 1.0;

	/* Assembly global system */
	uint32_t N_negative_jacobians = 0;
	for(uint32_t k=0; k < N_elements; k++){
		/* Get material properties */
		double density = vcn_fem_material_get_density(material);
		double E = vcn_fem_material_get_elasticity_module(material);

		if(elements_enabled != NULL){
			if(!elements_enabled[k]){
				E *= 1e-6;
				density *= 1e-6;
			}
		}

		double v = vcn_fem_material_get_poisson_module(material);
    
		/* Get constitutive matrix */
		double d11 = E/(1.0 - POW2(v));
		double d12 = v*d11;
    
		if(!enable_plane_stress){
			d11 = (E*(1.0-v))/((1.0 + v)*(1.0-2*v));
			d12 = (v*d11)/(1.0-v);
		}    
		double d22 = d11;
		double d33 = E/(2.0*(1.0+v));

		/* Allocate Cartesian derivatives for each Gauss Point */
		double* dNi_dx = (double*)malloc(elemtype->N_nodes*sizeof(double));
		double* dNi_dy = (double*)malloc(elemtype->N_nodes*sizeof(double));
		/* Compute constitutive matrix */
		double fx = 0.0;
		double fy = 0.0;
		if(enable_self_weight){
			fx = gravity[0] * vcn_fem_material_get_density(material);
			fy = gravity[1] * vcn_fem_material_get_density(material);
		}
    
		/* Integrate Ke and Fe using Gauss quadrature */
		memset(Ke, 0, 4 * elemtype->N_nodes * elemtype->N_nodes *
		       sizeof(double));
		if(M != NULL) 
			memset(Me, 0, 2 * elemtype->N_nodes * sizeof(double));
		memset(Fe, 0, 2 * elemtype->N_nodes * sizeof(double));
		for(uint32_t j=0; j < elemtype->N_Gauss_points; j++){      
			/* Get constitutive model */
			double dd11  = d11;
			double dd12  = d12;
			double dd22  = d22;
			double dd33  = d33;
      
			/* Compute Jacobian derivatives */
			double dx_dpsi = 0.0;
			double dy_dpsi = 0.0;
			double dx_deta = 0.0;
			double dy_deta = 0.0;
			for(uint32_t i=0; i < elemtype->N_nodes; i++){
				uint32_t inode = mesh->vertices_forming_triangles[k * elemtype->N_nodes + i];
				double xi = mesh->vertices[inode * 2];
				double yi = mesh->vertices[inode*2+1];
				dx_dpsi += 
					elemtype->dNi_dpsi[i](elemtype->psi[j], elemtype->eta[j]) * xi;
				dx_deta +=
					elemtype->dNi_deta[i](elemtype->psi[j], elemtype->eta[j]) * xi;
				dy_dpsi += 
					elemtype->dNi_dpsi[i](elemtype->psi[j], elemtype->eta[j]) * yi;
				dy_deta +=
					elemtype->dNi_deta[i](elemtype->psi[j], elemtype->eta[j]) * yi;
			}

			/* Compute Jacobian inverse and determinant */
			double detJ = dx_dpsi*dy_deta - dy_dpsi*dx_deta;
			/* Check if the element is distorted */
			if(detJ < 0) N_negative_jacobians += 1;
      
			double Jinv[4];
			Jinv[0] =  dy_deta/detJ;
			Jinv[1] = -dy_dpsi/detJ;
			Jinv[2] = -dx_deta/detJ;
			Jinv[3] =  dx_dpsi/detJ;
      
			/* Compute Shape functions derivatives in cartesian space */ 
			for(uint32_t i=0; i < elemtype->N_nodes; i++){
				dNi_dx[i] = 
					Jinv[0]*elemtype->dNi_dpsi[i](elemtype->psi[j], elemtype->eta[j]) + 
					Jinv[1]*elemtype->dNi_deta[i](elemtype->psi[j], elemtype->eta[j]);
				dNi_dy[i] = 
					Jinv[2]*elemtype->dNi_dpsi[i](elemtype->psi[j], elemtype->eta[j]) + 
					Jinv[3]*elemtype->dNi_deta[i](elemtype->psi[j], elemtype->eta[j]);
			}
			/* Compute elemental stiffness matrix
			 *       _            _         _           _        _  _
			 *      |dNi/dx        |       | d11 d12     |      |  |
			 * Bi = |       dNi/dy |   D = | d21 d22     |  K = |  | B'DB dx dy
			 *      |dNi/dy dNi/dx_|       |_        d33_|     _| _|

			 * B  = [B1 B2...Bi... Belemtype->N_nodes]
			 *
			 */
      
			for(uint32_t i1 = 0; i1 < elemtype->N_nodes; i1++){	
				for(uint32_t i2 = 0; i2 < elemtype->N_nodes; i2++){
					/*  Integrating elemental siffness matrix */
					Ke[(i1 * 2)*(2*elemtype->N_nodes) + (i2 * 2)] += 
						(dNi_dx[i1]*dNi_dx[i2]*dd11 + dNi_dy[i1]*dNi_dy[i2]*dd33) *
						detJ * thickness * elemtype->gp_weight[j];
					Ke[(i1 * 2)*(2*elemtype->N_nodes) + (i2*2+1)] +=
						(dNi_dx[i1]*dNi_dy[i2]*dd12 + dNi_dy[i1]*dNi_dx[i2]*dd33) *
						detJ * thickness * elemtype->gp_weight[j];
					Ke[(i1*2+1)*(2*elemtype->N_nodes) + (i2 * 2)] +=
						(dNi_dy[i1]*dNi_dx[i2]*dd12 + dNi_dx[i1]*dNi_dy[i2]*dd33) *
						detJ * thickness * elemtype->gp_weight[j];
					Ke[(i1*2+1)*(2*elemtype->N_nodes) + (i2*2+1)] +=
						(dNi_dy[i1]*dNi_dy[i2]*dd22 + dNi_dx[i1]*dNi_dx[i2]*dd33) *
						detJ * thickness * elemtype->gp_weight[j];
				}
				/* Calculate shape function of the i-th node at the j-th gauss point */
				double Ni_eval = 
					elemtype->Ni[i1](elemtype->psi[j], elemtype->eta[j]);
				/*  Integrating elemental mass matrix */
				if(M != NULL){
					/* OPPORTUNITY: Allocate just one for each component */
					Me[i1 * 2] += Ni_eval * Ni_eval * density *
						detJ * thickness * elemtype->gp_weight[j];
					Me[i1*2+1] += Ni_eval * Ni_eval * density * 
						detJ * thickness * elemtype->gp_weight[j];
				}
				/* Compute elemental forces vector */
				/* OPPORTUNITY: Allocate just one for each component */
				Fe[i1 * 2] += Ni_eval * fx * detJ * 
					thickness * elemtype->gp_weight[j];
				Fe[i1*2+1] += Ni_eval * fy * detJ * 
					thickness * elemtype->gp_weight[j];
			}
		}
		/* Add to global stiffness matrix */
		for(uint32_t i1 = 0; i1 < elemtype->N_nodes; i1++){
			for(uint32_t i2 = 0; i2 < elemtype->N_nodes; i2++){

				for(uint8_t j1=0; j1 < 2; j1++){
					for(uint8_t j2=0; j2 < 2; j2++){
						vcn_sparse_add
							(K, mesh->vertices_forming_triangles[k*3+i1]*2 + j1,
							 mesh->vertices_forming_triangles[k*3+i2]*2 + j2,
							 Ke[(i1*2+j1)*(2*elemtype->N_nodes) + (i2*2+j2)]);
					}
				}
			}
			/* Add to global mass matrix */
			if (M != NULL) {
				M[mesh->vertices_forming_triangles[k*3+i1] * 2] += Me[i1 * 2];
				M[mesh->vertices_forming_triangles[k*3+i1]*2+1] += Me[i1*2+1];
			}
			/* Add to global forces vector */
			F[mesh->vertices_forming_triangles[k*3+i1] * 2] += Fe[i1 * 2];
			F[mesh->vertices_forming_triangles[k*3+i1]*2+1] += Fe[i1*2+1];
		}
		free(dNi_dx);
		free(dNi_dy);
	}
	if (N_negative_jacobians > 0)
		/* OPPORTUNITY: Send to logfile */
		printf("ERROR: There are %i negative Jacobians.\n", N_negative_jacobians);

	/* Free elemental stiffness matrix and force vector */
	free(Ke);
	if (M != NULL) 
		free(Me);
	free(Fe);
}

void pipeline_set_boundary_conditions(vcn_sparse_t* K,
				      double* F, 
				      const vcn_bcond_t *const bmeshcond, 
				      double thickness,
				      double factor)
/* We assume that all the conditions are placed on the mesh nodes */
{
	/* OPPORTUNITY: Use libre_evaluator to eval time functions */
	/* Set Neuman conditions on vertices (Loads in the Solid Mechanics context) */
	for(uint32_t i=0; i < bmeshcond->N_Neuman_on_vtx; i++){
		for(uint32_t j=0; j < bmeshcond->N_dof; j++){
			if(!bmeshcond->Neuman_on_vtx_dof_mask[i * bmeshcond->N_dof + j])
				continue;
			F[bmeshcond->Neuman_on_vtx_idx[i] * bmeshcond->N_dof + j] += factor *
				bmeshcond->Neuman_on_vtx_val[i * bmeshcond->N_dof + j];
		}
	}

	/* Set Dirichlet conditions (Displacements in the Solid Mechanics context)  */
	for(uint32_t i=0; i < bmeshcond->N_Dirichlet_on_vtx; i++){
		for(uint32_t k=0; k < bmeshcond->N_dof; k++){
			if(!bmeshcond->Dirichlet_on_vtx_dof_mask[i * bmeshcond->N_dof + k])
				continue;
			uint32_t idx = bmeshcond->Dirichlet_on_vtx_idx[i] * bmeshcond->N_dof + k;
			double value = factor *
				bmeshcond->Dirichlet_on_vtx_val[i * bmeshcond->N_dof + k];
			vcn_sparse_set_Dirichlet_condition(K, F, idx, value);
		}
	}
}

void pipeline_compute_strain(double *strain,
			     const vcn_msh3trg_t *const mesh,
			     double *displacement,
			     const vcn_fem_elem_t *const elemtype,
			     bool enable_plane_stress,
			     const vcn_fem_material_t *const material)
{
	double *vertices = (double*) mesh->vertices;
	uint32_t N_elements = mesh->N_triangles;
	uint32_t *connectivity_mtx = (uint32_t*) mesh->vertices_forming_triangles;

	/* Initialize strains */
	memset(strain, 0, 3 * elemtype->N_Gauss_points * N_elements * sizeof(double));

	/* Iterate over elements to compute strain and stress at nodes */
	for(uint32_t k=0; k < N_elements; k++){
		double* dNi_dx = (double*)malloc(elemtype->N_nodes*sizeof(double));
		double* dNi_dy = (double*)malloc(elemtype->N_nodes*sizeof(double));

		/* Integrate domain */
		for(uint32_t j=0; j < elemtype->N_Gauss_points; j++){
			uint32_t idx = k * elemtype->N_Gauss_points + j;

			/* Compute Jacobian derivatives */
			double dx_dpsi = 0.0;
			double dy_dpsi = 0.0;
			double dx_deta = 0.0;
			double dy_deta = 0.0;

			for(uint32_t i=0; i < elemtype->N_nodes; i++){
				uint32_t inode = connectivity_mtx[k * elemtype->N_nodes + i];
				double xi = vertices[inode * 2];
				double yi = vertices[inode*2+1];
				dx_dpsi += 
					elemtype->dNi_dpsi[i](elemtype->psi[j], elemtype->eta[j]) * xi;
				dx_deta += 
					elemtype->dNi_deta[i](elemtype->psi[j], elemtype->eta[j]) * xi;
				dy_dpsi +=
					elemtype->dNi_dpsi[i](elemtype->psi[j], elemtype->eta[j]) * yi;
				dy_deta +=
					elemtype->dNi_deta[i](elemtype->psi[j], elemtype->eta[j]) * yi;
			}
      
			/* Compute Jacobian inverse and determinant */
			double detJ = dx_dpsi*dy_deta - dy_dpsi*dx_deta;
			double Jinv[4];
			Jinv[0] =  dy_deta/detJ;
			Jinv[1] = -dy_dpsi/detJ;
			Jinv[2] = -dx_deta/detJ;
			Jinv[3] =  dx_dpsi/detJ;
      
			/* Compute Shape functions derivatives in cartesian space */ 
			for(uint32_t i=0; i < elemtype->N_nodes; i++){
				dNi_dx[i] = 
					Jinv[0]*elemtype->dNi_dpsi[i](elemtype->psi[j], elemtype->eta[j]) + 
					Jinv[1]*elemtype->dNi_deta[i](elemtype->psi[j], elemtype->eta[j]);
				dNi_dy[i] = 
					Jinv[2]*elemtype->dNi_dpsi[i](elemtype->psi[j], elemtype->eta[j]) + 
					Jinv[3]*elemtype->dNi_deta[i](elemtype->psi[j], elemtype->eta[j]);
			}

			/* Compute Strain at Gauss Point */
			for (uint32_t i=0; i < elemtype->N_nodes; i++) {
				uint32_t inode = connectivity_mtx[k * elemtype->N_nodes + i];
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
				  const vcn_fem_elem_t *const elemtype)
{
	for(uint32_t i=0; i < N_elements; i++){
		for(uint32_t j=0; j < elemtype->N_Gauss_points; j++){
			uint32_t idx = i*elemtype->N_Gauss_points + j;
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
			 const vcn_fem_elem_t *const elemtype)
/* Compute elemental error based on strain field */
{
	uint32_t N_vertices = mesh->N_vertices;
	uint32_t N_elements = mesh->N_triangles;


	uint32_t N_gp = elemtype->N_Gauss_points;
	memset(error, 0, N_elements * N_gp * sizeof(double));

	/* Interpolate strain to nodes */
	double* strain_on_nodes = (double*)calloc(3 * N_vertices, sizeof(double));
	vcn_fem_interpolate_from_Gauss_points_to_nodes(mesh, elemtype,
						       3, strain, strain_on_nodes);

	uint32_t *connectivity_mtx = (uint32_t*) mesh->vertices_forming_triangles;

	/* Return values to Gauss Points to compute the error in the element */
	for(uint32_t i=0; i < N_elements; i++){
		error[i] = 0;
		for(uint32_t j=0; j < N_gp; j++){
			double strain_gp[3];
			memset(strain_gp, 0, 3*sizeof(double));
			for(uint32_t k=0; k < elemtype->N_nodes; k++){
				strain_gp[0] += elemtype->Ni[k](elemtype->psi[j], elemtype->eta[j]) * 
					strain_on_nodes[connectivity_mtx[i*3+k] * 3];
				strain_gp[1] += elemtype->Ni[k](elemtype->psi[j], elemtype->eta[j]) * 
					strain_on_nodes[connectivity_mtx[i*3+k]*3+1];
				strain_gp[2] += elemtype->Ni[k](elemtype->psi[j], elemtype->eta[j]) * 
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
