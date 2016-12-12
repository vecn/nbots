#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "nb/solver_bot.h"
#include "nb/memory_bot/allocate_mem.h"
#include "nb/geometric_bot.h"
#include "nb/pde_bot/material.h"
#include "nb/pde_bot/common_solid_mechanics/analysis2D.h"
#include "nb/pde_bot/common_solid_mechanics/formulas.h"
#include "nb/pde_bot/boundary_conditions/bcond.h"
#include "nb/pde_bot/boundary_conditions/bcond_iter.h"
#include "nb/pde_bot/finite_element/element.h"
#include "nb/pde_bot/finite_element/gaussp_to_nodes.h"
#include "nb/pde_bot/finite_element/solid_mechanics/static_damage2D.h"

#include "../utils.h"
#include "set_bconditions.h"
#include "pipeline.h"

#define MAX(a,b) (((a)>(b))?(a):(b))
#define POW2(a) ((a)*(a))

struct nb_fem_implicit_s{
	uint32_t N_steps;
	uint32_t N_max_iter;
	uint32_t N_max_iter_without_enhance;
	double residual_tolerance;
};
/*void damage_solver(nb_sparse_t *K, nb_sparse_t *L, nb_sparse_t *U,
                       uint32_t N_system_size, double *residual, double *displacement_increment, uint32_t solver_status,
                       bool enable_Cholesky_solver, uint32_t omp_parallel_threads);*/
int damage_solver(const nb_sparse_t *const A,
                   const double *const b, double* x);
static double tension_damage_r0(const nb_material_t *const mat);
static double tension_damage(const nb_material_t *const mat,
			     double *strain,
			     double *r_damage_prev,
			     double *r_damage,
			     double characteristic_length_of_fractured_domain,
			     nb_analysis2D_t analysis2D);
void nb_fem_compute_damage_stress_from_strain
			(uint32_t N_elements,
			 const nb_fem_elem_t *const elem,
			 const nb_material_t *const material,
			 nb_analysis2D_t analysis2D,
			 double* strain,
			 const bool* elements_enabled /* NULL to enable all */,
			 double* stress /* Output */,
			 double *damage_elem,
			 bool enable_computing_damage);
/*
SECOND OPTION
static double tension_truncated_damage_r0(const nb_material_t *const mat);
static double tension_truncated_damage
			(const nb_material_t *const mat,
			 double *strain,
			 double *r_damage_prev,
			 double* r_damage,
			 double characteristic_length_of_fractured_domain,
			 nb_analysis2D_t analysis2D);
*/
static void DMG_pipeline_assemble_system
		(nb_sparse_t* K, double* M, double *F,
		 const nb_mesh2D_t *const part,
		 const nb_fem_elem_t *const elem,
		 const nb_material_t *const material,
		 bool enable_self_weight,
		 double gravity[2],
		 nb_analysis2D_t analysis2D,
		 nb_analysis2D_params *params2D,
		 bool enable_computing_damage,
		 double* damage_elem,
		 bool* elements_enabled /* NULL to enable all */);

static void DMG_pipeline_compute_strain
			(double *strain,
			 const nb_mesh2D_t *const part,
			 double *displacement,
			 const nb_fem_elem_t *const elem,
			 bool enable_computing_damage,
			 nb_analysis2D_t analysis2D,
			 const nb_material_t *const material,
			 double *damage,
			 double *r_dmg_prev,
			 double *r_dmg);

static double get_clfd(const nb_mesh2D_t *part, uint32_t id_elem);

nb_fem_implicit_t* nb_fem_implicit_create(void)
{
	return nb_allocate_zero_mem(sizeof(nb_fem_implicit_t));
}

void nb_fem_implicit_destroy(nb_fem_implicit_t* isparams){
	nb_free_mem(isparams);
}

void nb_fem_implicit_set_N_steps(nb_fem_implicit_t* isparams,
				  uint32_t N_steps){
	isparams->N_steps = N_steps;
}

void nb_fem_implicit_set_N_max_iter(nb_fem_implicit_t* isparams,
				     uint32_t N_max_iter){
	isparams->N_max_iter = N_max_iter;
}

void nb_fem_implicit_set_N_max_iter_without_enhance
(nb_fem_implicit_t* isparams, uint32_t N_max_iter){
	isparams->N_max_iter_without_enhance = N_max_iter;
}

void nb_fem_implicit_set_residual_tolerance
(nb_fem_implicit_t* isparams, double penergy_tol)
{
	isparams->residual_tolerance = penergy_tol;
}

uint32_t nb_fem_implicit_get_N_steps(nb_fem_implicit_t* isparams){
	return isparams->N_steps;
}

uint32_t nb_fem_implicit_get_N_max_iter(nb_fem_implicit_t* isparams){
	return isparams->N_max_iter;
}

uint32_t nb_fem_implicit_get_N_max_iter_without_enhance
(nb_fem_implicit_t* isparams){
	return isparams->N_max_iter_without_enhance;
}

double nb_fem_implicit_get_residual_tolerance
(nb_fem_implicit_t* isparams)
{
	return isparams->residual_tolerance;
}

static inline double tension_damage_r0(const nb_material_t *const mat)
{
	return nb_material_get_traction_limit_stress(mat) /
		sqrt(nb_material_get_elasticity_module(mat));
}

static double tension_damage(const nb_material_t *const mat,
			     double *strain,
			     double *r_damage_prev,
			     double *r_damage,
			     double characteristic_length_of_fractured_domain,
			     nb_analysis2D_t analysis2D)
{
	double D[4];
	nb_pde_get_constitutive_matrix(D, mat, analysis2D);
	double E = nb_material_get_elasticity_module(mat);
	double Gf = nb_material_get_fracture_energy(mat);
	double ft = nb_material_get_traction_limit_stress(mat);

	/* Compute effective stress */
	double effective_stress[3];
	effective_stress[0] = D[0] * strain[0] + D[1] * strain[1];
	effective_stress[1] = D[1] * strain[0] + D[2] * strain[1];
	effective_stress[2] = D[3] * strain[2];
	/* Compute principal stress using Mohr's circle */
	double sigma_avg = 0.5 * (effective_stress[0] + effective_stress[1]);
	double R = sqrt(0.25 *
			POW2(effective_stress[0] - effective_stress[1]) +
			POW2(effective_stress[2]));

	double positive_stress1 = MAX(0, sigma_avg + R);
	double positive_stress2 = MAX(0, sigma_avg - R);

	/* Compute inverse of D */
	double detD2x2 = D[0]*D[2] - D[1]*D[1];
	double id11 =  D[2]/detD2x2;
	double id12 = -D[1]/detD2x2;
	double id22 =  D[0]/detD2x2;
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
		(Gf*E)/(characteristic_length_of_fractured_domain*POW2(ft));
	double A = 1.0 / (div - 0.5);
	double G = 1.0 - (r0/r_damage[0])*exp(A*(1.0-(r_damage[0]/r0)));
	return G;
	}

/*
SECOND OPTION
static inline double tension_truncated_damage_r0
			(const nb_material_t *const mat)
{
	return nb_material_get_traction_limit_stress(mat);
}

static double tension_truncated_damage
			(const nb_material_t *const mat,
			 double *strain,
			 double *r_damage_prev,
			 double* r_damage,
			 double characteristic_length_of_fractured_domain,
			 nb_analysis2D_t analysis2D)
{
	double D[4];
	nb_pde_get_constitutive_matrix(D, mat, analysis2D);
	double E = nb_material_get_elasticity_module(mat);
	double Gf = nb_material_get_fracture_energy(mat);
	double ft = nb_material_get_traction_limit_stress(mat);

// Compute effective stress
	double effective_stress[3];
	effective_stress[0] = D[0] * strain[0] + D[1] * strain[1];
	effective_stress[1] = D[1] * strain[0] + D[2] * strain[1];
	effective_stress[2] = D[3] * strain[2] * 0.5;

// Compute principal stress using Mohr's circle
	double sigma_avg = 0.5 * (effective_stress[0] + effective_stress[1]);
	double R = sqrt(0.25 *
			POW2(effective_stress[0] - effective_stress[1]) +
			POW2(effective_stress[2]));
	double main_stress1 = sigma_avg + R;
	double main_stress2 = sigma_avg - R;

// Compute tau
	double tau = MAX(0, main_stress1) +
		MAX(0, main_stress2);

// Compute and return damage
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
*/

uint8_t nb_fem_compute_2D_Damage_Solid_Mechanics
			(const nb_mesh2D_t *const part,
			 const nb_fem_elem_t *const elem,
			 const nb_material_t *const material,
			 const nb_bcond_t *const bcond,
			 bool enable_self_weight,
			 double gravity[2],
			 bool enable_Cholesky_solver,
			 nb_analysis2D_t analysis2D,
			 nb_analysis2D_params *params2D,
			 nb_fem_implicit_t* params,
			 const char* logfile,
			 double *damage,
			 double *displacement,
			 double *strain,
			 double *stress)
/* Quasistatic formulation */
{
    uint8_t status = 1;
	uint64_t N_nod = nb_mesh2D_get_N_nodes(part);
	uint64_t N_elem = nb_mesh2D_get_N_elems(part);

	uint32_t omp_parallel_threads = 1;

	FILE *log = fopen(logfile, "a");
	fprintf(log, "FEM: Damage Model\n");
	fclose(log);

	uint64_t N_system_size = 2*N_nod;

	uint8_t N_gp = nb_fem_elem_get_N_gpoints(elem);

	/*******************************************************************/
	/*********************** > ?????? **********************************/
	/*******************************************************************/
	memset(displacement, 0, 2 * N_nod * sizeof(*displacement));
	memset(strain, 0, 3 * N_elem * sizeof(*strain));
	memset(damage, 0, N_elem * N_gp * sizeof(double));
	memset(stress, 0, 3 * N_elem * sizeof(*stress));
	double* r_dmg = nb_allocate_mem(N_gp * N_elem * sizeof(double));

	/* Initialize r parameter used for damage calculation */
	for (uint32_t i = 0; i < N_gp * N_elem; i++)
		r_dmg[i] = tension_damage_r0(material);

	/*******************************************************************/
	/****************** > Allocate system ******************************/
	/*******************************************************************/
	/* Allocate global Stiffness Matrices */
	nb_graph_t *graph = malloc(nb_graph_get_memsize());
	nb_graph_init(graph);
	nb_mesh2D_load_graph(part, graph, NB_NODES_LINKED_BY_ELEMS);
	nb_sparse_t *K = nb_sparse_create(graph, NULL, 2);
	nb_graph_finish(graph);
	nb_sparse_t *L = NULL;
	nb_sparse_t *U = NULL;
	/* Allocate the triangular matrices L and U using symbolic Cholesky */
	//nb_sparse_alloc_LU(K, &L, &U);


	/* Allocate force vectors and displacement increment */
	double* F = nb_allocate_zero_mem(N_system_size*sizeof(double));
	double* P = nb_allocate_zero_mem(N_system_size*sizeof(double));
	double* residual = nb_allocate_zero_mem(N_system_size*sizeof(double));
	double* displacement_increment = nb_allocate_zero_mem(N_system_size*sizeof(double));

	/* Allocate damage parameter 'r' */
	double* r_dmg_prev = nb_allocate_mem(N_gp * N_elem * sizeof(double));
    bool react_known_node = false;
    double max_y_force = 0;
	/*******************************************************************/
	/******************* > Start simulation of N steps *****************/
	/*******************************************************************/
	for (uint32_t n = 0; n < nb_fem_implicit_get_N_steps(params); n++) {
        printf("Step %d\n", n);
		log = fopen(logfile, "a");
		fprintf(log, "  [ Load step %i]\n", n + 1);
		fclose(log);
		memcpy(r_dmg_prev, r_dmg, N_gp * N_elem * sizeof(double));

		/***********************************************************/
		/*************** > Implicit integration ********************/
		/***********************************************************/
		/* Implicit integration scheme */
		if(!enable_Cholesky_solver)
			memset(displacement_increment, 0, N_system_size * sizeof(double));

		double residual_norm = 1;
		uint32_t residual_iter = 0;
		uint32_t residual_iter_without_enhance = 0;
		double residual_best;
        uint32_t solver_status = 0; /* Default initialization (Must be initialized) */
        memset(F, 0, N_system_size);
		while(1) {
			/**************************************************/
			/********* > Assemble system **********************/
			/**************************************************/
			DMG_pipeline_assemble_system(K, NULL, F,
						     part,
						     elem,
						     material,
						     enable_self_weight,
						     gravity,
						     analysis2D,
						     params2D,
						     true, /* Enable computing damage */
						     damage,
						     NULL);

			/*****************************************/
			/****** > Set boundary conditions ********/
			/*****************************************/
			double condition_factor =
				(n + 1.0)/(double) nb_fem_implicit_get_N_steps(params);

			/* Set Boundary Conditions */
			nb_fem_set_bconditions(part, K, F, bcond, condition_factor);

			/*******************************************/
			/******* > Verify residual *****************/
			/*******************************************/
			/* Compute P increment */
			nb_sparse_multiply_vector(K, displacement, P, omp_parallel_threads);

			/* Compute residual norm */
			residual_norm = 0;
			for(int i = 0; i < N_system_size; i++){
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
			solver_status = damage_solver(K, residual, displacement_increment);
			/*damage_solver(K, L, U, N_system_size, residual, displacement_increment, solver_status,
                                 enable_Cholesky_solver, omp_parallel_threads);*/
            if(solver_status != 0){
                goto END_OF_PROCESS;
            }
			/* Write logfile */
			log = fopen(logfile, "a");
			fprintf(log, "    [%i: Residual: %e]\n", residual_iter, residual_norm);
			fclose(log);

			/* Update displacements, velocities and accelerations */
			for(int i=0; i < N_system_size; i++){
				displacement[i] += displacement_increment[i];
			}

			/**********************************************/
			/***** > Compute Strain and Damage      *******/
			/**********************************************/
			DMG_pipeline_compute_strain(strain, part, displacement, elem, true,
						    analysis2D,
						    material, damage,
						    r_dmg_prev, r_dmg);

			/* Increase iterator */
			residual_iter ++;
		}
	}

    nb_fem_compute_damage_stress_from_strain(N_elem, elem, material, analysis2D, strain,
                                             NULL, stress, damage, true);

	/*****************************************************************/
	/******************** > Free memory ******************************/
	/*****************************************************************/
	status = 0;

	END_OF_PROCESS:
	nb_sparse_destroy(K);
	if(enable_Cholesky_solver){
		nb_sparse_destroy(L);
		nb_sparse_destroy(U);
	}

	nb_free_mem(F);
	nb_free_mem(displacement_increment);
	nb_free_mem(residual);
	nb_free_mem(P);

	nb_free_mem(r_dmg_prev);
	nb_free_mem(r_dmg);

	return status;
}

static void DMG_pipeline_assemble_system
		(nb_sparse_t* K, double* M, double *F,
		 const nb_mesh2D_t *const part,
		 const nb_fem_elem_t *const elem,
		 const nb_material_t *const material,
		 bool enable_self_weight,
		 double gravity[2],
		 nb_analysis2D_t analysis2D,
		 nb_analysis2D_params *params2D,
		 bool enable_computing_damage,
		 double* damage_elem,
		 bool* elements_enabled /* NULL to enable all */)
{
	uint32_t N_elems = nb_mesh2D_get_N_elems(part);
	uint32_t N_nod = nb_mesh2D_get_N_nodes(part);

	nb_sparse_reset(K);
	if (NULL != M)
		memset(M, 0, nb_sparse_get_size(K) * sizeof(double));
	memset(F, 0, 2*N_nod*sizeof(double));

	/* Allocate elemental Stiffness Matrix and Force Vector */
	uint8_t N_nodes = nb_fem_elem_get_N_nodes(elem);
	double* Ke = nb_allocate_mem(4 * POW2(N_nodes) * sizeof(double));
	double* Me = NULL;
	if(M != NULL)
		Me = nb_allocate_mem(2 * N_nodes * sizeof(double));
	double* Fe = nb_allocate_mem(2 * N_nodes * sizeof(double));

	/* Assembly global system */
	uint32_t N_negative_jacobians = 0;
	for (uint32_t k = 0; k < N_elems; k++) {
		double D[4] = {1e-6, 1e-6, 1e-6, 1e-6};
		double density = 1e-6;
		if (pipeline_elem_is_enabled(elements_enabled, k)) {
            density = nb_material_get_density(material);
			nb_pde_get_constitutive_matrix(D, material,
						       analysis2D);
		}

		/* Allocate Cartesian derivatives for each Gauss Point */
		double* dNi_dx = nb_allocate_mem(N_nodes * sizeof(double));
		double* dNi_dy = nb_allocate_mem(N_nodes * sizeof(double));

		/* Compute constitutive matrix */
		double fx = 0.0;
		double fy = 0.0;
		if(enable_self_weight){
			fx = gravity[0] * density;
			fy = gravity[1] * density;
		}

		/* Integrate Ke and Fe using Gauss quadrature */
		memset(Ke, 0, 4 * POW2(N_nodes) * sizeof(double));
		if(M != NULL)
			memset(Me, 0, 2 * N_nodes * sizeof(double));
		memset(Fe, 0, 2 * N_nodes * sizeof(double));

		uint8_t N_gp = nb_fem_elem_get_N_gpoints(elem);
		for (uint32_t j = 0; j < N_gp; j++) {
			/* Get constitutive model */
			double Dr[4];
			memcpy(Dr, D, 4 * sizeof(double));
			if (enable_computing_damage) {
				Dr[0] *= (1.0 - damage_elem[k * N_gp + j]);
				Dr[1] *= (1.0 - damage_elem[k * N_gp + j]);
				Dr[2] *= (1.0 - damage_elem[k * N_gp + j]);
				Dr[3] *= (1.0 - damage_elem[k * N_gp + j]);
			}

			/* Compute Jacobian derivatives */
			double Jinv[4];
			double detJ = nb_fem_get_jacobian(elem, k, part, j, Jinv);

			if (0 > detJ)
				goto EXIT;

			nb_fem_get_derivatives(elem, j, Jinv, dNi_dx, dNi_dy);

			double thickness = params2D->thickness;
			pipeline_sum_gauss_point(elem, j, Dr, density, thickness,
						 detJ, dNi_dx, dNi_dy, fx, fy,
						 Ke, Me, Fe);
		}

		pipeline_add_to_global_system(elem, k, part, Ke, Me, Fe,
					      K, M, F);

		nb_free_mem(dNi_dx);
		nb_free_mem(dNi_dy);
	}
 EXIT:
	/* Free elemental stiffness matrix and force vector */
	nb_free_mem(Ke);
	if (M != NULL)
		nb_free_mem(Me);
	nb_free_mem(Fe);
}

static void DMG_pipeline_compute_strain
			(double *strain,
			 const nb_mesh2D_t *const part,
			 double *displacement,
			 const nb_fem_elem_t *const elem,
			 bool enable_computing_damage,
			 nb_analysis2D_t analysis2D,
			 const nb_material_t *const material,
			 double *damage,
			 double *r_dmg_prev,
			 double *r_dmg)
{
	uint32_t N_elems = nb_mesh2D_get_N_elems(part);

	/* Initialize strains */
	uint8_t N_gp = nb_fem_elem_get_N_gpoints(elem);
	memset(strain, 0, 3 * N_gp * N_elems * sizeof(double));
    double max_damage = 0.0;
	/* Iterate over elements to compute strain, stress and damage at nodes */
	for (uint32_t k = 0 ; k < N_elems; k++) {
		uint8_t N_nodes = nb_fem_elem_get_N_nodes(elem);
		double* dNi_dx = nb_allocate_mem(N_nodes * sizeof(double));
		double* dNi_dy = nb_allocate_mem(N_nodes * sizeof(double));

		/* Integrate domain */
		uint8_t N_gp = nb_fem_elem_get_N_gpoints(elem);
		for (uint32_t j = 0; j < N_gp; j++) {
			double Jinv[4];
			double detJ = nb_fem_get_jacobian(elem, k, part, j, Jinv);
			if (nb_fem_elem_is_distorted(detJ))
				goto EXIT;

			nb_fem_get_derivatives(elem, j, Jinv, dNi_dx, dNi_dy);

			uint32_t idx = k * N_gp + j;
			/* Compute Strain at Gauss Point */
			for (uint32_t i = 0; i < N_nodes; i++) {
				uint32_t inode = nb_mesh2D_elem_get_adj(part, k, i);
				strain[idx * 3] += dNi_dx[i] * displacement[inode * 2];
				strain[idx*3+1] += dNi_dy[i] * displacement[inode*2+1];
				strain[idx*3+2] += (dNi_dy[i] * displacement[inode * 2] +
						    dNi_dx[i] * displacement[inode*2+1]);
			}
			/* Compute damage */
			if (enable_computing_damage) {
				double clfd = get_clfd(part, k);

				damage[idx] = tension_damage(material,
							     &(strain[idx*3]),
							     &(r_dmg_prev[idx]),
							     &(r_dmg[idx]),
							     clfd,
							     analysis2D);
                if(damage[idx] > max_damage) /*TEMPORAL*/
                    max_damage = damage[idx];
			}
		}
		nb_free_mem(dNi_dx);
		nb_free_mem(dNi_dy);
	}
EXIT:
	return;
}

static double get_clfd(const nb_mesh2D_t *part, uint32_t id_elem)
/* characteristic_length_of_fractured_domain */
{
	uint32_t n1 = nb_mesh2D_elem_get_adj(part, id_elem, 0);
	uint32_t n2 = nb_mesh2D_elem_get_adj(part, id_elem, 1);
	uint32_t n3 = nb_mesh2D_elem_get_adj(part, id_elem, 2);
	double x1 = nb_mesh2D_node_get_x(part, n1);
	double y1 = nb_mesh2D_node_get_y(part, n1);
	double x2 = nb_mesh2D_node_get_x(part, n2);
	double y2 = nb_mesh2D_node_get_y(part, n2);
	double x3 = nb_mesh2D_node_get_x(part, n3);
	double y3 = nb_mesh2D_node_get_y(part, n3);
	return sqrt((x2 - x1) * (y3 - y1) -
		    (y2 - y1) * (x3 - x1));
}
int damage_solver(const nb_sparse_t *const A,
		  const double *const b, double* x)
{
	uint32_t N = nb_sparse_get_size(A);
	memset(x, 0, N * sizeof(*x));
    double err;
        uint32_t iter;
	int status = nb_sparse_solve_CG_precond_Jacobi(A, b, x, N,
							1e-8, &iter,
							&err, 1);
	int out;
	if (0 == status || 1 == status)
		out = 0;
	else {
        printf("Tolerance not reached in CG Jacobi.\n");
        printf("Error: %lf\n",err);
        printf("Iterations: %d\n", iter);
		out = 1; // Tolerance not reached in CG Jacobi
    }
	return out;
}

//void damage_solver(const nb_sparse_t *K, nb_sparse_t *L, nb_sparse_t *U, uint32_t N_system_size, double *residual, double *displacement_increment,
                   //uint32_t solver_status, bool enable_Cholesky_solver, uint32_t omp_parallel_threads) {
   //if(enable_Cholesky_solver){
        /* Decompose matrix */
        //solver_status =
          //  nb_sparse_decompose_Cholesky(K, L, U, omp_parallel_threads);

        //if(solver_status != 0)
           // nb_sparse_decompose_LU(K, L, U, omp_parallel_threads);

            /* Solve system */
        //nb_sparse_solve_LU(L, U, residual, displacement_increment);
    //}

   // if(!enable_Cholesky_solver || solver_status != 0){
        //if(solver_status != 0){
            /*log = fopen(logfile, "a");
            fprintf(log, "    [Cholesky Fails (The stiffness matrix is not positive definite).]\n");
            fprintf(log, "    [Solving with Conjugate Gradient (Tolerance: 1e-8)...]\n");
            fclose(log);*/
      //  }
        //double error = 0.0; /* Default initialization (Must be initialized) */
        //uint32_t iters = 0;     /* Default initialization (Must be initialized) */
        //memset(displacement_increment, 0, N_system_size*sizeof(double));
        //uint64_t N = nb_sparse_get_size(K);
        //solver_status =
          //  nb_sparse_solve_CG_precond_Jacobi(K, residual, displacement_increment, N, 1e-8,
            //                    NULL, NULL, 1);
        //if(solver_status != 0){
          //  printf("Precond Jacobi solver failed!\n");
            //solver_status = nb_sparse_solve_Gauss_Seidel(K, residual, displacement_increment,
              //                          60 * N, 1e-8, NULL, NULL, 1);
            //if(solver_status != 0){
            //printf("Gauss seidel solver failed!\n");
                /*log = fopen(logfile, "a");
                fprintf(log, "    [Conjugate Gradient (p. Jacobi) and Gauss-Seidel fails.]\n");
                fprintf(log, "      [Gauss-Seidel iterations: %u]\n", iters);
                fprintf(log, "      [Gauss-Seidel error: %e ]\n", error);
                fclose(log);*/
        //    }
        //}
    //}
//}

void nb_fem_compute_damage_stress_from_strain
			(uint32_t N_elements,
			 const nb_fem_elem_t *const elem,
			 const nb_material_t *const material,
			 nb_analysis2D_t analysis2D,
			 double* strain,
			 const bool* elements_enabled /* NULL to enable all */,
			 double* stress /* Output */,
			 double *damage_elem,
			 bool enable_computing_damage)
{
	/* Compute stress from element strain */
	uint32_t omp_parallel_threads = 1;
	uint8_t N_gp = nb_fem_elem_get_N_gpoints(elem);
#pragma omp parallel for num_threads(omp_parallel_threads) schedule(guided)
	for (uint32_t i = 0; i < N_elements; i++) {
		double D[4] = {1e-6, 1e-6, 1e-6, 1e-6};
		if (pipeline_elem_is_enabled(elements_enabled, i))
			nb_pde_get_constitutive_matrix(D, material,
						       analysis2D);
		for (uint32_t j = 0; j < N_gp; j++) {
			/* Get constitutive model */
			double Dr[4];
			memcpy(Dr, D, 4 * sizeof(double));
			if (enable_computing_damage) {
				Dr[0] *= (1.0 - damage_elem[i * N_gp + j]);
				Dr[1] *= (1.0 - damage_elem[i * N_gp + j]);
				Dr[2] *= (1.0 - damage_elem[i * N_gp + j]);
				Dr[3] *= (1.0 - damage_elem[i * N_gp + j]);
			}
			uint32_t id = i * N_gp + j;
			stress[id * 3] = strain[id * 3] * Dr[0] +
				strain[id*3+1] * Dr[1];
			stress[id*3+1] = strain[id * 3] * Dr[1] +
				strain[id*3+1] * Dr[2];
			stress[id*3+2] = strain[id*3+2] * Dr[3];
        }
    }
}
