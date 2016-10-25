#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "nb/memory_bot.h"
#include "nb/solver_bot/sparse/sparse.h"
#include "nb/solver_bot/sparse/solvers/cg_precond_jacobi.h"

#include "../sparse_struct.h"

int nb_sparse_solve_CG_precond_Jacobi
(const nb_sparse_t *const A, 
 const double *const b,
 double *_x,                /* Out */
 uint32_t max_iter, double tolerance,
 uint32_t* niter_performed,     /* Out (NULL if not required) */
 double* tolerance_reached, /* Out (NULL if not required) */
 uint32_t omp_parallel_threads)
/* Return the num of iterations */
{
	/* Solve Ax = b with Conjugate Gradient preconditioned with jacobi */
	uint32_t vec_size = A->N * sizeof(double);
	char *memblock = nb_allocate_zero_mem(5 * vec_size);	
	double* g = (void*) memblock;
	double* p = (void*) (memblock + vec_size);
	double* q = (void*) (memblock + 2 * vec_size);
	double* w = (void*) (memblock + 3 * vec_size);
	double* Aii = (void*) (memblock + 4 * vec_size);
	double dot_gg = 0;

#pragma omp parallel for reduction(+:dot_gg) num_threads(omp_parallel_threads) schedule(guided)
	for (uint32_t i = 0; i < A->N; i++) {
		double sum = 0;
		for(uint32_t j=0; j< A->rows_size[i]; j++)
			sum += A->rows_values[i][j] * _x[A->rows_index[i][j]];
		g[i] = sum - b[i];
		Aii[i] = nb_sparse_get(A,i,i);
		q[i] = g[i]/Aii[i];
		p[i] = -q[i];
		dot_gg += g[i]*g[i];
	}
	uint32_t k = 0;
	while (dot_gg > tolerance * tolerance && k < max_iter) {
		double dot_pw = 0;
		double dot_gq = 0;
		dot_gg = 0;

#pragma omp parallel for reduction(+:dot_pw, dot_gg, dot_gq) num_threads(omp_parallel_threads)
		for (uint32_t i = 0; i < A->N; i++) {
			w[i] = 0;
			for (uint32_t j = 0; j <  A->rows_size[i]; j++)
				w[i] += A->rows_values[i][j] * p[A->rows_index[i][j]];
			dot_pw += p[i]*w[i];
			dot_gg += g[i]*g[i];
			dot_gq += g[i]*q[i];
		}
		double alphak = dot_gq/dot_pw;
		double dot_gkqk = 0;
		
#pragma omp parallel for reduction(+:dot_gkqk) num_threads(omp_parallel_threads)
		for (uint32_t i=0; i< A->N; i++) {
			_x[i] += alphak*p[i];
			g[i] += alphak*w[i];
			q[i] = g[i]/Aii[i];
			dot_gkqk += g[i]*q[i];
		}

		double betak = dot_gkqk/dot_gq;
		
#pragma omp parallel for num_threads(omp_parallel_threads)
		for (uint32_t i=0; i< A->N; i++)
			p[i] = -q[i]+betak*p[i];
		k++;
	}
	/* Free memory */
	nb_free_mem(memblock);

	if (NULL != niter_performed)
		niter_performed[0]= k;

	if (NULL != tolerance_reached)
		*tolerance_reached = sqrt(dot_gg);

	if (dot_gg > tolerance*tolerance)
		return 1;

	return 0;
}
