struct vcn_fem_implicit_s{
	uint32_t N_steps;
	uint32_t N_max_iter;
	uint32_t N_max_iter_without_enhance;
	double residual_tolerance;
};
vcn_fem_implicit_t* vcn_fem_implicit_create(){
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

double tension_damage_r0(const vcn_fem_material_t *const mat){
	return vcn_fem_material_get_traction_limit_stress(mat) /
		sqrt(vcn_fem_material_get_elasticity_module(mat));
}

double tension_damage
(const vcn_fem_material_t *const mat, 
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
	double d11 = E/(1.0 - vcn_math_pow2(v));
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
			vcn_math_pow2(effective_stress[0] - effective_stress[1]) +
			vcn_math_pow2(effective_stress[2]));
	double positive_stress1 = vcn_math_maxd(0, sigma_avg + R);
	double positive_stress2 = vcn_math_maxd(0, sigma_avg - R);

	/* Compute inverse of D */
	double detD2x2 = d11*d22 - d12*d12;
	double id11 =  d22/detD2x2;
	double id12 = -d12/detD2x2;
	double id22 =  d11/detD2x2;
	/* Compute stress ^T strain */
	double sTs =
		id11 * vcn_math_pow2(positive_stress1) + 
		2 * id12 * positive_stress1 * positive_stress2 +
		id22 * vcn_math_pow2(positive_stress2);

	/* Compute tau */
	double tau = sqrt(sTs);

	/* Compute and return damage */
	double r0 = mat->r0_damage(mat);
	r_damage[0] = vcn_math_maxd(r_damage_prev[0], tau);
	double div =  
		(Gf/characteristic_length_of_fractured_domain)*(E/vcn_math_pow2(ft));
	double A = 1.0 / (div - 0.5);
	double G = 1.0 - (r0/r_damage[0])*exp(A*(1.0-(r_damage[0]/r0)));
	return G;
}

double tension_truncated_damage_r0(const vcn_fem_material_t *const mat){
	return vcn_fem_material_get_traction_limit_stress(mat);
}

double tension_truncated_damage
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
	double d11 = E/(1.0 - vcn_math_pow2(v));
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
			vcn_math_pow2(effective_stress[0] - effective_stress[1]) +
			vcn_math_pow2(effective_stress[2]));
	double main_stress1 = sigma_avg + R;
	double main_stress2 = sigma_avg - R;

	/* Compute tau */
	double tau = vcn_math_maxd(0, main_stress1) + 
		vcn_math_maxd(0, main_stress2);

	/* Compute and return damage */
	double r0 = mat->r0_damage(mat);
	double beta = 2.0;
	r_damage[0] = vcn_math_maxd(r_damage_prev[0], tau);
	double pi = (beta * characteristic_length_of_fractured_domain *
		     vcn_math_pow2(ft))
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
 vcn_fem_output_t* output,
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

	/****************************************************************************/
	/*********************** > Allocate output **********************************/
	/****************************************************************************/
	double* displacement = (double*)calloc(2 * N_vertices, sizeof(double));
	double* strain = (double*)calloc(3 * N_elements * N_gp, sizeof(double));
	double* damage = (double*)calloc(N_elements * N_gp, sizeof(double));
	double* r_dmg = (double*)malloc(N_gp * N_elements * sizeof(double));

	/* Initialize r parameter used for damage calculation */
	for (uint32_t i = 0; i < N_gp * N_elements; i++)
		r_dmg[i] = material->r0_damage(material);

	/****************************************************************************/
	/*********************** > Allocate system **********************************/
	/****************************************************************************/
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
  
	/************************************************************************/
	/************************ > Start simulation of N steps *****************/
	/************************************************************************/
	for (uint32_t n = 0; n < vcn_fem_implicit_get_N_steps(params); n++) {
		log = fopen(logfile, "a");
		fprintf(log, "  [ Load step %i]\n", n + 1);
		fclose(log);
		memcpy(r_dmg_prev, r_dmg, N_gp * N_elements * sizeof(double));

		/***********************************************************************/
		/******************** > Implicit integration ***************************/
		/***********************************************************************/
		/* Implicit integration scheme */
		if(!enable_Cholesky_solver)
			memset(du, 0, N_system_size * sizeof(double));

		double residual_norm = 1;
		uint32_t residual_iter = 0;
		uint32_t residual_iter_without_enhance = 0;
		double residual_best;
		while(1) {
			/**********************************************************************/
			/******************** > Assemble system *******************************/
			/**********************************************************************/
			pipeline_assemble_system(K, NULL, F,
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

			/**********************************************************************/
			/******************** > Set boundary conditions ***********************/
			/**********************************************************************/
			double condition_factor =
				(n + 1.0)/(double) vcn_fem_implicit_get_N_steps(params);

			/* Set Boundary Conditions */
			pipeline_set_boundary_conditions(K, F, bmeshcond, thickness, condition_factor);

			/**********************************************************************/
			/******************** > Verify residual *******************************/
			/**********************************************************************/
			/* Compute P increment */
			vcn_sparse_multiply_vector(K, displacement, P, omp_parallel_threads);

			/* Compute residual norm */
			residual_norm = 0;
			for(uint32_t i=0; i < N_system_size; i++){
				residual[i] = F[i] - P[i];
				residual_norm += vcn_math_pow2(residual[i]);
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

			/**********************************************************************/
			/************** > Solve system (to compute displacements) *************/
			/**********************************************************************/
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
            

			/**********************************************************************/
			/**************** > Compute Strain and Damage         *****************/
			/**********************************************************************/
			pipeline_compute_strain(strain, mesh, displacement, elemtype, true,
						enable_plane_stress, material, damage,
						r_dmg_prev, r_dmg);

			/* Increase iterator */
			residual_iter ++;
		}
	}
	/*****************************************************************************/
	/***************************** > Free memory *********************************/
	/*****************************************************************************/
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
