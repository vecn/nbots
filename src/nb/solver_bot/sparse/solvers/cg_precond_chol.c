#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "nb/math_bot.h"
#include "nb/memory_bot.h"
#include "nb/container_bot.h"
#include "nb/solver_bot/sparse/sparse.h"
#include "nb/solver_bot/sparse/solvers/lu.h"
#include "nb/solver_bot/sparse/solvers/cholesky.h"
#include "nb/solver_bot/sparse/solvers/cg_precond_chol.h"

#include "../sparse_struct.h"
#include "cholesky_symbolic.h"

int nb_sparse_solve_CG_precond_Cholesky(const nb_sparse_t *const A,
					const double *const b, 
					double *_x,                /* Out */
					uint32_t ktrunc,
					uint32_t max_iter, double tolerance,
					uint32_t* niter_performed,     /* Out (NULL if not required) */
					double* tolerance_reached, /* Out (NULL if not required) */
					uint32_t omp_parallel_threads)
/* Return the num of iterations */
{
	/* Solve Ax = b with Conjugate Gradient preconditioned with 
	 * Cholesky truncated.
	 */
	double* g = nb_allocate_zero_mem(A->N * sizeof(double));
	double* p = nb_allocate_zero_mem(A->N * sizeof(double));
	double* q = nb_allocate_zero_mem(A->N * sizeof(double));
	double* w = nb_allocate_zero_mem(A->N * sizeof(double));

	double dot_gg = 0;
  
#pragma omp parallel for reduction(+:dot_gg) num_threads(omp_parallel_threads) schedule(guided)
	for(uint32_t i=0; i< A->N; i++){
		double sum = 0;
		for(uint32_t j=0; j< A->rows_size[i]; j++)
			sum += A->rows_values[i][j] * _x[A->rows_index[i][j]];
		g[i] = sum - b[i];
		dot_gg += g[i]*g[i];
	}
	/* Solve H Ht q = g for q */
	nb_sparse_t *H = nb_sparse_allocate(A->N);
	nb_sparse_t *Ht = nb_sparse_allocate(A->N);
	nb_sparse_cholesky_symbolic(A, H, Ht, ktrunc);

	nb_sparse_decompose_Cholesky(A,H,Ht, omp_parallel_threads);
	nb_sparse_solve_LU(H,Ht,g,q);

#pragma omp parallel for num_threads(omp_parallel_threads)
	for(uint32_t i=0; i< A->N; i++)
		p[i] = -q[i];

	uint32_t k = 0;
	while(dot_gg > tolerance*tolerance && k < max_iter){
		double dot_pw = 0;
		double dot_gq = 0;
		dot_gg = 0;
#pragma omp parallel for reduction(+:dot_pw, dot_gg, dot_gq) num_threads(omp_parallel_threads)
		for(uint32_t i=0; i < A->N; i++){
			w[i] = 0;
			for(uint32_t j=0; j< A->rows_size[i]; j++)
				w[i] += A->rows_values[i][j] * p[A->rows_index[i][j]];
			dot_pw += p[i]*w[i];
			dot_gg += g[i]*g[i];
			dot_gq += g[i]*q[i];
		}
		double alphak = dot_gq/dot_pw;
#pragma omp parallel for num_threads(omp_parallel_threads)
		for(uint32_t i=0; i < A->N; i++){
			_x[i] += alphak*p[i];
			g[i] += alphak*w[i];
		}
		/* Solve H Ht q = g for q */
		nb_sparse_solve_LU(H,Ht,g,q);

		double dot_gkqk = 0;
#pragma omp parallel for reduction(+:dot_gkqk) num_threads(omp_parallel_threads)
		for(uint32_t i=0; i< A->N; i++)
			dot_gkqk += g[i]*q[i];

		double betak = dot_gkqk/dot_gq;
#pragma omp parallel for num_threads(omp_parallel_threads)
		for(uint32_t i=0; i< A->N; i++)
			p[i] = -q[i] + betak*p[i];
		k++;
	}
	/* Free memory */
	nb_free_mem(g);
	nb_free_mem(p);
	nb_free_mem(q);
	nb_free_mem(w);
	nb_sparse_destroy(H);
	nb_sparse_destroy(Ht);

	if(niter_performed != NULL) niter_performed[0]= k;

	if(tolerance_reached != NULL) *tolerance_reached = sqrt(dot_gg);

	if(dot_gg > tolerance*tolerance)
		return 1;

	return 0;
}
