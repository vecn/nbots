#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "nb/memory_bot.h"
#include "nb/solver_bot/sparse/sparse.h"
#include "nb/solver_bot/sparse/solvers/conjugate_gradient.h"

#include "../sparse_struct.h"

int nb_sparse_solve_conjugate_gradient
(const nb_sparse_t *const A, 
 const double *const b, 
 double *_x,                /* Out */
 uint32_t max_iter, double tolerance,
 uint32_t* niter_performed,     /* Out (NULL if not required) */
 double* tolerance_reached, /* Out (NULL if not required) */
 uint32_t omp_parallel_threads)
/* Return the num of iterations */
{
	/* Solve Ax = b with Conjugate Gradient method */
	double *g = nb_allocate_zero_mem(A->N * sizeof(double));
	double *p = nb_allocate_zero_mem(A->N * sizeof(double));
	double *w = nb_allocate_zero_mem(A->N * sizeof(double));
	double dot_gg = 0;

#pragma omp parallel for reduction(+:dot_gg) num_threads(omp_parallel_threads) schedule(guided)
	for(uint32_t i=0; i< A->N; i++){
		double sum = 0;
		for(uint32_t j=0; j< A->rows_size[i]; j++)
			sum += A->rows_values[i][j] * _x[A->rows_index[i][j]];
		g[i] = sum - b[i];
		p[i] = -g[i];
		dot_gg += g[i]*g[i];
	}
	uint32_t k = 0;
	while(dot_gg > tolerance*tolerance && k < max_iter){
		double dot_pw = 0;
		dot_gg = 0;
#pragma omp parallel for reduction(+:dot_pw, dot_gg) num_threads(omp_parallel_threads) schedule(guided)
		for(uint32_t i = 0; i< A->N; i++){
			w[i] = 0;
			for(uint32_t j = 0; j< A->rows_size[i]; j++)
				w[i] += A->rows_values[i][j] * p[A->rows_index[i][j]];
			dot_pw += p[i]*w[i];
			dot_gg += g[i]*g[i];
		}
		double alphak = dot_gg/dot_pw;
		double dot_gkgk = 0;
#pragma omp parallel for reduction(+:dot_gkgk) num_threads(omp_parallel_threads) schedule(guided)
		for(uint32_t i=0; i< A->N; i++){
			_x[i] += alphak*p[i];
			g[i] += alphak*w[i];
			dot_gkgk += g[i]*g[i];
		}
		double betak = dot_gkgk/dot_gg;
#pragma omp parallel for num_threads(omp_parallel_threads)
		for(uint32_t i=0; i< A->N; i++)
			p[i] = -g[i] + betak * p[i];
		k++;
	}
	/* Free memory */
	nb_free_mem(g);
	nb_free_mem(p);
	nb_free_mem(w);
  
	if(niter_performed != NULL) niter_performed[0]= k;

	if(tolerance_reached != NULL) *tolerance_reached = sqrt(dot_gg);

	if(dot_gg > tolerance*tolerance)
		return 1;

	return 0;
}
