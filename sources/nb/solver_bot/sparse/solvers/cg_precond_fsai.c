#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "nb/math_bot.h"
#include "nb/memory_bot.h"
#include "nb/container_bot.h"
#include "nb/solver_bot.h"

#include "../sparse_struct.h"

#define POW2(a) ((a)*(a))

int nb_sparse_solve_CG_precond_fsai
(const nb_sparse_t *const A,
 const double *const b,
 double *_x,                /* Out */
 double threshold,
 uint32_t max_iter,	double tolerance,
 uint32_t* niter_performed,     /* Out (NULL if not required) */
 double* tolerance_reached, /* Out (NULL if not required) */
 uint32_t omp_parallel_threads)
/* Return the num of iterations */
{
	/* Conjugate gradient preconditioned with "Factorized sparse 
	 * approximated inverse" 
	 */
	double *D = nb_allocate_zero_mem(A->N * sizeof(double));
	double *siD = nb_allocate_zero_mem(A->N * sizeof(double));

	nb_sparse_t* G  = nb_sparse_allocate(A->N);
	nb_sparse_t* Gt = nb_sparse_allocate(A->N);
	/* Generate D diagonal matrix as
	 *           
	 *     Dii = |Aii|,   if |Aii| > 0
	 *            1       otherwise
	 */

#pragma omp parallel for num_threads(omp_parallel_threads)
	for(uint32_t i=0; i < A->N; i++){
		D[i] = nb_sparse_get(A,i,i);
		if(D[i] == 0)
			D[i] = 1;
		/* Compute D^(-1/2) */
		siD[i] = 1/sqrt(D[i]);
	}
  
	/* Generate structure of G lower triangular matrix 
	 *
	 *    G = 1 ,  if (i == j) or (|[D^(-1/2) A D^(-1/2)]_ij| > threshold)
	 *        0 ,  otherwise
	 *
	 */
	for(uint32_t i=0; i < A->N; i++){
		uint32_t isize = 0;
		uint32_t isizet = 0;

#pragma omp parallel for reduction(+:isize,isizet) num_threads(omp_parallel_threads)
		for(uint32_t q=0; q< A->rows_size[i]; q++){
			uint32_t j = A->rows_index[i][q];
			if(i == j || 
			   fabs(siD[i] * A->rows_values[i][q] * siD[j])
			   > threshold){
				if(i > j)
					isize++;
				else if(i < j)
					isizet++;
				else{
					isize++;
					isizet++;
				}
			}
		}
		G->rows_size[i] = isize;
		G->rows_index[i] = nb_allocate_zero_mem(isize *
							sizeof(uint32_t));
		G->rows_values[i] = nb_allocate_zero_mem(isize *
							 sizeof(double));

		Gt->rows_size[i] = isizet;
		Gt->rows_index[i] = nb_allocate_zero_mem(isizet *
							 sizeof(uint32_t));
		Gt->rows_values[i] = nb_allocate_zero_mem(isizet *
							  sizeof(double));
	}

#pragma omp parallel for num_threads(omp_parallel_threads)
	for(uint32_t i=0; i < A->N; i++){
		/* Compute values of ~G */
		double* subA =
			nb_allocate_zero_mem(POW2(G->rows_size[i]) *
					     sizeof(double));
		/* The data of vector g is not allocated, is a pointer to each row of ~G */
		double* subg = G->rows_values[i];
		double *delta = nb_allocate_zero_mem(G->rows_size[i] *
						     sizeof(double));
		uint32_t k = 0;
		for(uint32_t q = 0; q < A->rows_size[i]; q++){
			uint32_t j = A->rows_index[i][q];
			if(i == j || 
			   fabs(siD[i] * A->rows_values[i][q] * siD[j]) > threshold){
				if(i >= j){
					G->rows_index[i][k] = j;
					for(uint32_t l=0; l<k; l++){
						subA[k*G->rows_size[i] + l] = 
							nb_sparse_get(A,j,G->rows_index[i][l]);
						subA[l*G->rows_size[i] + k] = 
							nb_sparse_get(A,G->rows_index[i][l],j);
					}
					subA[k*G->rows_size[i] + k] = nb_sparse_get(A,j,j);
					if(i == j)
						delta[k] = 1;
					else
						delta[k] = 0;
					k++;
				}
			}
		}
		double* L = nb_allocate_zero_mem(POW2(G->rows_size[i]) *
						 sizeof(double));
		nb_matrix_cholesky_decomposition(subA, L, G->rows_size[i]);
		nb_matrix_cholesky_solve(L, delta, subg, G->rows_size[i]);
		/* Finally do G = [~G]*D   */
		for(uint32_t q=0; q < G->rows_size[i]; q++)
			G->rows_values[i][q] *= D[G->rows_index[i][q]];

		/* Free memory */
		nb_free_mem(subA);
		nb_free_mem(L);
		nb_free_mem(delta);
	}
	/* Store G transposed */
	nb_sparse_get_transpose(G,Gt);

	/* Free memory */
	nb_free_mem(D);
	nb_free_mem(siD);

	/* Solve Ax = b with Conjugate Gradient method */
	double* r = nb_allocate_zero_mem(A->N * sizeof(double));
	double* p = nb_allocate_zero_mem(A->N * sizeof(double));
	double* w = nb_allocate_zero_mem(A->N * sizeof(double));
	double* Gr = nb_allocate_zero_mem(A->N * sizeof(double));
	double* Mr = nb_allocate_zero_mem(A->N * sizeof(double));

	double dot_rr = 0;

#pragma omp parallel for reduction(+:dot_rr) schedule(guided) num_threads(omp_parallel_threads)
	for(uint32_t i=0; i< A->N; i++){
		double Ax_i = 0;
		for(uint32_t j=0; j< A->rows_size[i]; j++)
			Ax_i += A->rows_values[i][j] * _x[A->rows_index[i][j]];
		r[i] = b[i] - Ax_i;
		dot_rr += r[i]*r[i];
	}
	/* Compute Gr */

#pragma omp parallel for num_threads(omp_parallel_threads)
	for(uint32_t i=0; i< A->N; i++){
		Gr[i] = 0;
		for(uint32_t j=0; j< G->rows_size[i]; j++)
			Gr[i] += G->rows_values[i][j] * r[G->rows_index[i][j]];
	}
	/* Compute Mr <- G'(Gr) */

#pragma omp parallel for num_threads(omp_parallel_threads)
	for(uint32_t i=0; i< A->N; i++){
		Mr[i] = 0;
		for(uint32_t j=0; j< Gt->rows_size[i]; j++)
			Mr[i] += Gt->rows_values[i][j] * Gr[Gt->rows_index[i][j]];
		p[i] = Mr[i];
	}
	uint32_t k = 0;
	/* Start iterations */
	while(dot_rr > tolerance*tolerance && k < max_iter){
		double dot_pw = 0;
		double dot_rMr = 0;

#pragma omp parallel for reduction(+:dot_pw, dot_rMr) num_threads(omp_parallel_threads) schedule(guided)
		for(uint32_t i=0; i< A->N; i++){
			w[i] = 0;
			for(uint32_t j=0; j< A->rows_size[i]; j++)
				w[i] += A->rows_values[i][j] * p[A->rows_index[i][j]];
			dot_pw += p[i]*w[i];
			dot_rMr += r[i]*Mr[i];
		}
		double alphak = dot_rMr/dot_pw;
		dot_rr = 0;

#pragma omp parallel for reduction(+:dot_rr) num_threads(omp_parallel_threads) schedule(guided)
		for(uint32_t i=0; i< A->N; i++){
			_x[i] += alphak*p[i];
			r[i] -= alphak*w[i];
			dot_rr += r[i]*r[i];
		}
		/* Compute Gr */

#pragma omp parallel for num_threads(omp_parallel_threads)
		for(uint32_t i=0; i< A->N; i++){
			Gr[i] = 0;
			for(uint32_t j=0; j< G->rows_size[i]; j++)
				Gr[i] += G->rows_values[i][j] * r[G->rows_index[i][j]];
		}
		/* Compute Mr <- G'(Gr) */
		double dot_rkMrk = 0;

#pragma omp parallel for reduction(+:dot_rkMrk) num_threads(omp_parallel_threads)
		for(uint32_t i=0; i< A->N; i++){
			Mr[i] = 0;
			for(uint32_t j=0; j< Gt->rows_size[i]; j++)
				Mr[i] += Gt->rows_values[i][j] * Gr[Gt->rows_index[i][j]];
			dot_rkMrk += r[i]*Mr[i];
		}
		double betak = dot_rkMrk/dot_rMr;

#pragma omp parallel for num_threads(omp_parallel_threads)
		for(uint32_t i=0; i< A->N; i++)
			p[i] = Mr[i] + betak*p[i];
		k++;
	}
	/* Free memory */
	nb_sparse_destroy(G);
	nb_sparse_destroy(Gt);
	nb_free_mem(r);
	nb_free_mem(p);
	nb_free_mem(w);
	nb_free_mem(Gr);
	nb_free_mem(Mr);

	if(niter_performed != NULL) niter_performed[0]= k;

	if(tolerance_reached != NULL) *tolerance_reached = sqrt(dot_rr);

	if(dot_rr > tolerance*tolerance)
		return 1;

	return 0;
}
