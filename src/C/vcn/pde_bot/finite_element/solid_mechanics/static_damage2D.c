#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "vcn/eigen_bot.h"
#include "vcn/pde_bot/material.h"
#include "vcn/pde_bot/boundary_conditions.h"
#include "vcn/pde_bot/finite_element/element.h"
#include "vcn/pde_bot/finite_element/gaussp_to_nodes.h"
#include "vcn/pde_bot/finite_element/solid_mechanics/static_damage2D.h"
  
#include "../element_struct.h"
#include "pipeline.h"

#define MAX(a,b) (((a)>(b))?(a):(b))
#define POW2(a) ((a)*(a))

struct vcn_fem_implicit_s{
	uint32_t N_steps;
	uint32_t N_max_iter;
	uint32_t N_max_iter_without_enhance;
	double residual_tolerance;
};

static double tension_damage_r0(const vcn_fem_material_t *const mat);
static double tension_damage(const vcn_fem_material_t *const mat, 
			     double *strain,
			     double *r_damage_prev,
			     double *r_damage,
			     double characteristic_length_of_fractured_domain,
			     bool enable_plane_stress);
static double tension_truncated_damage_r0(const vcn_fem_material_t *const mat);
static double tension_truncated_damage
			(const vcn_fem_material_t *const mat, 
			 double *strain,
			 double *r_damage_prev,
			 double* r_damage,
			 double characteristic_length_of_fractured_domain,
			 bool enable_plane_stress);
static void DMG_pipeline_assemble_system
		(vcn_sparse_t* K, double* M, double *F,
		 const vcn_msh3trg_t *const mesh,
		 const vcn_fem_elem_t *const elemtype,
		 const vcn_fem_material_t *const material,
		 bool enable_self_weight,
		 double gravity[2],
		 bool enable_plane_stress,
		 double thickness,
		 bool enable_computing_damage,
		 double* damage_elem,
		 bool* elements_enabled /* NULL to enable all */);

static void DMG_pipeline_compute_strain
			(double *strain,
			 const vcn_msh3trg_t *const mesh,
			 double *displacement,
			 const vcn_fem_elem_t *const elemtype,
			 bool enable_computing_damage,
			 bool enable_plane_stress,
			 const vcn_fem_material_t *const material,
			 double *damage,
			 double *r_dmg_prev,
			 double *r_dmg);

vcn_fem_implicit_t* vcn_fem_implicit_create(void)
{
	return (vcn_fem_implicit_t*)
		calloc(1, sizeof(vcn_fem_implicit_t));
}

void vcn_fem_implicit_destroy(vcn_fem_implicit_t* isparams){
	free(isparams);
}

void vcn_fem_implicit_set_N_steps(vcn_fem_implicit_t* isparams,
				  uint32_t N_steps){
	isparams->N_steps = N_steps;
}

void vcn_fem_implicit_set_N_max_iter(vcn_fem_implicit_t* isparams,
				     uint32_t N_max_iter){
	isparams->N_max_iter = N_max_iter;
}

void vcn_fem_implicit_set_N_max_iter_without_enhance
(vcn_fem_implicit_t* isparams, uint32_t N_max_iter){
	isparams->N_max_iter_without_enhance = N_max_iter;
}

void vcn_fem_implicit_set_residual_tolerance
(vcn_fem_implicit_t* isparams, double penergy_tol)
{
	isparams->residual_tolerance = penergy_tol;
}

uint32_t vcn_fem_implicit_get_N_steps(vcn_fem_implicit_t* isparams){
	return isparams->N_steps;
}

uint32_t vcn_fem_implicit_get_N_max_iter(vcn_fem_implicit_t* isparams){
	return isparams->N_max_iter;
}

uint32_t vcn_fem_implicit_get_N_max_iter_without_enhance
(vcn_fem_implicit_t* isparams){
	return isparams->N_max_iter_without_enhance;
}

double vcn_fem_implicit_get_residual_tolerance
(vcn_fem_implicit_t* isparams)
{
	return isparams->residual_tolerance;
}

static inline double tension_damage_r0(const vcn_fem_material_t *const mat)
{
	return vcn_fem_material_get_traction_limit_stress(mat) /
		sqrt(vcn_fem_material_get_elasticity_module(mat));
}

static double tension_damage(const vcn_fem_material_t *const mat, 
			     double *strain,
			     double *r_damage_prev,
			     double *r_damage,
			     double characteristic_length_of_fractured_domain,
			     bool enable_plane_stress)
{
	/* Material properties */
	double E = vcn_fem_material_get_elasticity_module(mat);
	double v = vcn_fem_material_get_poisson_module(mat);
	double Gf = vcn_fem_material_get_fracture_energy(mat);
	double ft = vcn_fem_material_get_traction_limit_stress(mat);

	/* Get constitutive matrix */
	double d11 = E/(1.0 - POW2(v));
	double d12 = v*d11;
    
	if(!enable_plane_stress){
		d11 = (E*(1.0-v))/((1.0 + v)*(1.0-2*v));
		d12 = (v*d11)/(1.0-v);
	}    
	double d22 = d11;
	double d33 = E/(2.0*(1.0+v));

	/* Compute effective stress */
	double effective_stress[3];
	effective_stress[0] = d11 * strain[0] + d12 * strain[1];
	effective_stress[1] = d12 * strain[0] + d22 * strain[1];
	effective_stress[2] = d33 * strain[2];

	/* Compute principal stress using Mohr's circle */
	double sigma_avg = 0.5 * (effective_stress[0] + effective_stress[1]);
	double R = sqrt(0.25 * 
			POW2(effective_stress[0] - effective_stress[1]) +
			POW2(effective_stress[2]));
	double positive_stress1 = MAX(0, sigma_avg + R);
	double positive_stress2 = MAX(0, sigma_avg - R);

	/* Compute inverse of D */
	double detD2x2 = d11*d22 - d12*d12;
	double id11 =  d22/detD2x2;
	double id12 = -d12/detD2x2;
	double id22 =  d11/detD2x2;
	/* Compute stress ^T strain */
	double sTs =
		id11 * POW2(positive_stress1) + 
		2 * id12 * positive_stress1 * positive_stress2 +
		id22 * POW2(positive_stress2);

	/* Compute tau */
	double tau = sqrt(sTs);

	/* Compute and return damage */
	double r0 = tension_damage_r0(mat);
	r_damage[0] = MAX(r_damage_prev[0], tau);
	double div =  
		(Gf/characteristic_length_of_fractured_domain)*(E/POW2(ft));
	double A = 1.0 / (div - 0.5);
	double G = 1.0 - (r0/r_damage[0])*exp(A*(1.0-(r_damage[0]/r0)));
	return G;
}

static inline double tension_truncated_damage_r0
			(const vcn_fem_material_t *const mat)
{
	return vcn_fem_material_get_traction_limit_stress(mat);
}

static double tension_truncated_damage
			(const vcn_fem_material_t *const mat, 
			 double *strain,
			 double *r_damage_prev,
			 double* r_damage,
			 double characteristic_length_of_fractured_domain,
			 bool enable_plane_stress)
{
	/* Material properties */
	double E = vcn_fem_material_get_elasticity_module(mat);
	double v = vcn_fem_material_get_poisson_module(mat);
	double Gf = vcn_fem_material_get_fracture_energy(mat);
	double ft = vcn_fem_material_get_traction_limit_stress(mat);

	/* Get constitutive matrix */
	double d11 = E/(1.0 - POW2(v));
	double d12 = v*d11;
    
	if(!enable_plane_stress){
		d11 = (E*(1.0-v))/((1.0 + v)*(1.0-2*v));
		d12 = (v*d11)/(1.0-v);
	}    
	double d22 = d11;
	double d33 = E/(2.0*(1.0+v));

	/* Compute effective stress */
	double effective_stress[3];
	effective_stress[0] = d11 * strain[0] + d12 * strain[1];
	effective_stress[1] = d12 * strain[0] + d22 * strain[1];
	effective_stress[2] = d33 * strain[2] * 0.5;

	/* Compute principal stress using Mohr's circle */
	double sigma_avg = 0.5 * (effective_stress[0] + effective_stress[1]);
	double R = sqrt(0.25 * 
			POW2(effective_stress[0] - effective_stress[1]) +
			POW2(effective_stress[2]));
	double main_stress1 = sigma_avg + R;
	double main_stress2 = sigma_avg - R;

	/* Compute tau */
	double tau = MAX(0, main_stress1) + 
		MAX(0, main_stress2);

	/* Compute and return damage */
	double r0 = tension_damage_r0(mat);
	double beta = 2.0;
	r_damage[0] = MAX(r_damage_prev[0], tau);
	double pi = (beta * characteristic_length_of_fractured_domain *
		     POW2(ft))
		/(2*E*Gf);
	if(pi > 1.0) pi = 1.0;
	double Hs = pi/(1.0 - pi);
	double G = 1.0 - (r0/r_damage[0])*exp(2.0*Hs*(1.0-r_damage[0]/r0));
	return G;
}

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
			 bool restore_computation,
			 const char* logfile)
/* Quasistatic formulation */
{
	uint32_t N_vertices = mesh->N_vertices;
	uint32_t N_elements = mesh->N_triangles;
	double* vertices = mesh->vertices;
	uint32_t* connectivity_mtx = mesh->vertices_forming_triangles;

	uint32_t omp_parallel_threads = 1;

	FILE *log = fopen(logfile, "a");
	fprintf(log, "FEM: Damage Model\n");
	fclose(log);
  
	uint32_t N_system_size = N_vertices * 2;

	uint32_t N_gp = elemtype->N_Gauss_points;

	/*******************************************************************/
	/*********************** > ?????? **********************************/
	/*******************************************************************/
	double* displacement = (double*)calloc(2 * N_vertices, sizeof(double));
	double* strain = (double*)calloc(3 * N_elements * N_gp, sizeof(double));
	double* damage = (double*)calloc(N_elements * N_gp, sizeof(double));
	double* r_dmg = (double*)malloc(N_gp * N_elements * sizeof(double));

	/* Initialize r parameter used for damage calculation */
	for (uint32_t i = 0; i < N_gp * N_elements; i++)
		r_dmg[i] = tension_damage_r0(material);

	/*******************************************************************/
	/****************** > Allocate system ******************************/
	/*******************************************************************/
	/* Allocate global Stiffness Matrices */
	vcn_graph_t *graph = vcn_msh3trg_create_vtx_graph(mesh);
	vcn_sparse_t* K = vcn_sparse_create(graph, NULL, 2);
	vcn_sparse_t *L = NULL;
	vcn_sparse_t *U = NULL;
	/* Allocate the triangular matrices L and U using symbolic Cholesky */
	vcn_sparse_alloc_LU(K, &L, &U);

  
	/* Allocate force vectors and displacement increment */
	double* F = (double*)calloc(N_system_size, sizeof(double));
	double* P = (double*)calloc(N_system_size, sizeof(double));
	double* residual = (double*)calloc(N_system_size, sizeof(double));
	double* du = (double*)calloc(N_system_size, sizeof(double));

	/* Allocate damage parameter 'r' */
	double* r_dmg_prev = (double*)malloc(N_gp * N_elements * sizeof(double));
  
	/*******************************************************************/
	/******************* > Start simulation of N steps *****************/
	/*******************************************************************/
	for (uint32_t n = 0; n < vcn_fem_implicit_get_N_steps(params); n++) {
		log = fopen(logfile, "a");
		fprintf(log, "  [ Load step %i]\n", n + 1);
		fclose(log);
		memcpy(r_dmg_prev, r_dmg, N_gp * N_elements * sizeof(double));

		/***********************************************************/
		/*************** > Implicit integration ********************/
		/***********************************************************/
		/* Implicit integration scheme */
		if(!enable_Cholesky_solver)
			memset(du, 0, N_system_size * sizeof(double));

		double residual_norm = 1;
		uint32_t residual_iter = 0;
		uint32_t residual_iter_without_enhance = 0;
		double residual_best;
		while(1) {
			/**************************************************/
			/********* > Assemble system **********************/
			/**************************************************/
			DMG_pipeline_assemble_system(K, NULL, F,
						 mesh,
						 elemtype,
						 material,
						 enable_self_weight,
						 gravity,
						 enable_plane_stress,
						 thickness,
						 true, /* Enable computing damage */
						 damage,
						 NULL);

			/*****************************************/
			/****** > Set boundary conditions ********/
			/*****************************************/
			double condition_factor =
				(n + 1.0)/(double) vcn_fem_implicit_get_N_steps(params);

			/* Set Boundary Conditions */
			pipeline_set_boundary_conditions(K, F, bmeshcond, thickness, condition_factor);

			/*******************************************/
			/******* > Verify residual *****************/
			/*******************************************/
			/* Compute P increment */
			vcn_sparse_multiply_vector(K, displacement, P, omp_parallel_threads);

			/* Compute residual norm */
			residual_norm = 0;
			for(uint32_t i=0; i < N_system_size; i++){
				residual[i] = F[i] - P[i];
				residual_norm += POW2(residual[i]);
			}
      
			residual_norm = sqrt(residual_norm);

			/* Check if residual is minimized */
			if(residual_iter == 0){
				residual_best = residual_norm;
			}else{
				if(residual_norm < residual_best){
					residual_best = residual_norm;
					residual_iter_without_enhance = 0;
				}else
					residual_iter_without_enhance += 1;
			}

			/* Check stop criteria */
			if((residual_norm < params->residual_tolerance
			    && residual_iter > 0) || residual_iter >= params->N_max_iter ||
			   residual_iter_without_enhance >= params->N_max_iter_without_enhance){
				/* Write logfile */
				log = fopen(logfile, "a");
				fprintf(log, "    [%i: Residual: %e ]\n",
					residual_iter, residual_norm);
				fclose(log);
				break;
			}

			/***********************************************/
			/** > Solve system (to compute displacements) **/
			/***********************************************/
			char solver_status = 0; /* Default initialization (Must be initialized) */
			if(enable_Cholesky_solver){
				/* Decompose matrix */
				solver_status = 
					vcn_sparse_decompose_Cholesky(K, L, U, omp_parallel_threads);

				if(solver_status != 0)
					vcn_sparse_decompose_LU(K, L, U, omp_parallel_threads);
  
				/* Solve system */
				vcn_sparse_solve_LU(L, U, residual, du);
			}
			if(!enable_Cholesky_solver || solver_status != 0){
				if(solver_status != 0){
					log = fopen(logfile, "a");
					fprintf(log, "    [Cholesky Fails (The stiffness matrix is not positive definite).]\n");
					fprintf(log, "    [Solving with Conjugate Gradient (Tolerance: 1e-8)...]\n");
					fclose(log);
				}
				double error = 0.0; /* Default initialization (Must be initialized) */
				uint32_t iters = 0;     /* Default initialization (Must be initialized) */
				solver_status = 
					vcn_sparse_solve_CG_precond_Jacobi(K, residual, du,
									   vcn_sparse_get_size(K),
									   1e-8,
									   &iters, &error,
									   omp_parallel_threads);
				if(solver_status != 0){
					solver_status = vcn_sparse_solve_Gauss_Seidel(K, residual, du,
										      60 * vcn_sparse_get_size(K),
										      1e-8,
										      &iters, &error,
										      omp_parallel_threads);
					if(solver_status != 0){
						log = fopen(logfile, "a");
						fprintf(log, "    [Conjugate Gradient (p. Jacobi) and Gauss-Seidel fails.]\n");
						fprintf(log, "      [Gauss-Seidel iterations: %u]\n", iters);
						fprintf(log, "      [Gauss-Seidel error: %e ]\n", error);
						fclose(log);
					}
				}
			}
			/* Write logfile */
			log = fopen(logfile, "a");
			fprintf(log, "    [%i: Residual: %e]\n", residual_iter, residual_norm);
			fclose(log);

			/* Update displacements, velocities and accelerations */
			for(uint32_t i=0; i < N_system_size; i++)
				displacement[i] += du[i];
            

			/**********************************************/
			/***** > Compute Strain and Damage      *******/
			/**********************************************/
			DMG_pipeline_compute_strain(strain, mesh, displacement, elemtype, true,
						enable_plane_stress, material, damage,
						r_dmg_prev, r_dmg);

			/* Increase iterator */
			residual_iter ++;
		}
	}
	/*****************************************************************/
	/******************** > Free memory ******************************/
	/*****************************************************************/
	free(displacement);
	free(strain);
	free(damage);

	vcn_sparse_destroy(K);
	if(enable_Cholesky_solver){
		vcn_sparse_destroy(L);
		vcn_sparse_destroy(U);
	}

	free(F);
	free(du);
	free(residual);
	free(P); 
  
	free(r_dmg_prev);
	free(r_dmg);

	if(restore_computation){
		free(vertices);
		free(connectivity_mtx);
	}
}

static void DMG_pipeline_assemble_system
		(vcn_sparse_t* K, double* M, double *F,
		 const vcn_msh3trg_t *const mesh,
		 const vcn_fem_elem_t *const elemtype,
		 const vcn_fem_material_t *const material,
		 bool enable_self_weight,
		 double gravity[2],
		 bool enable_plane_stress,
		 double thickness,
		 bool enable_computing_damage,
		 double* damage_elem,
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
      
			if(enable_computing_damage){
				dd11 *= (1.0 - damage_elem[k*elemtype->N_Gauss_points + j]);
				dd12 *= (1.0 - damage_elem[k*elemtype->N_Gauss_points + j]);
				dd22 *= (1.0 - damage_elem[k*elemtype->N_Gauss_points + j]);
				dd33 *= (1.0 - damage_elem[k*elemtype->N_Gauss_points + j]);
			}

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

static void DMG_pipeline_compute_strain
			(double *strain,
			 const vcn_msh3trg_t *const mesh,
			 double *displacement,
			 const vcn_fem_elem_t *const elemtype,
			 bool enable_computing_damage,
			 bool enable_plane_stress,
			 const vcn_fem_material_t *const material,
			 double *damage,
			 double *r_dmg_prev,
			 double *r_dmg)
{
	double *vertices = (double*) mesh->vertices;
	uint32_t N_elements = mesh->N_triangles;
	uint32_t *connectivity_mtx = (uint32_t*) mesh->vertices_forming_triangles;

	/* Initialize strains */
	memset(strain, 0, 3 * elemtype->N_Gauss_points * N_elements * sizeof(double));

	/* Iterate over elements to compute strain, stress and damage at nodes */
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
			for(uint32_t i=0; i < elemtype->N_nodes; i++){
				uint32_t inode = connectivity_mtx[k * elemtype->N_nodes + i];
				strain[idx * 3] += dNi_dx[i] * displacement[inode * 2];
				strain[idx*3+1] += dNi_dy[i] * displacement[inode*2+1];
				strain[idx*3+2] += (dNi_dy[i] * displacement[inode * 2] +
						    dNi_dx[i] * displacement[inode*2+1]);
			}
      
			/* Compute damage */
			if(enable_computing_damage){
				uint32_t n1 = connectivity_mtx[k * 3];
				uint32_t n2 = connectivity_mtx[k*3+1];
				uint32_t n3 = connectivity_mtx[k*3+2];
				double characteristic_length_of_fractured_domain =
					sqrt(((vertices[n2 * 2]-vertices[n1 * 2]) *
					      (vertices[n3*2+1]-vertices[n1*2+1]) -
					      (vertices[n2*2+1]-vertices[n1*2+1]) *
					      (vertices[n3 * 2]-vertices[n1 * 2])));
      
				damage[idx] = 
					tension_damage(material, 
							 &(strain[idx*3]),
							 &(r_dmg_prev[idx]),
							 &(r_dmg[idx]),
							 characteristic_length_of_fractured_domain,
							 enable_plane_stress);
			}
		}
		free(dNi_dx);
		free(dNi_dy);
	}
}
