#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "nb/memory_bot.h"
#include "nb/solver_bot.h"

#include "../sparse_struct.h"

#define POW2(a) ((a)*(a))

int nb_sparse_solve_Gauss_Seidel
			(const nb_sparse_t *const A, 
			 const double *const b,
			 double *_x,                /* Out */
			 uint32_t max_iter, double tolerance,
			 uint32_t* niter_performed,     /* Out (NULL if not required) */
			 double* tolerance_reached, /* Out (NULL if not required) */
			 uint32_t omp_parallel_threads)
{
	/* Allocate RHS vector */
	double* c = (double*)nb_allocate_mem(nb_sparse_get_size(A)*sizeof(double));
	/* Start iterations */
	register uint32_t k = 0;
	double error = 1e10;
	while(k < max_iter) {
		/* Generate RHS vector */
		error = 0;
#pragma omp parallel for reduction(+:error) num_threads(omp_parallel_threads)
		for(uint32_t i=0; i < nb_sparse_get_size(A); i++){
			c[i] = b[i];
			double error_i = b[i];
			for(uint32_t j=0; j < A->rows_size[i]; j++){
				error_i -= A->rows_values[i][j] * _x[A->rows_index[i][j]];
				/* Use only the Upper triangular matrix to create the RHS */
				if(A->rows_index[i][j] > i)
					c[i] -= A->rows_values[i][j] * _x[A->rows_index[i][j]];
			}
			error += POW2(error_i);
		}
		/* Check convergence */
		if(error <= tolerance*tolerance)
			break;

		/* Solve by forward substitution */
		nb_sparse_forward_solve(A, c, _x);

		/* Increase iterations */
		k++;
	}
	/* Free memory */
	nb_free_mem(c);

	/* Successful exit */
	if(niter_performed != NULL) *niter_performed = k;

	if(tolerance_reached != NULL)
		*tolerance_reached = sqrt(error);
  
	if(error > tolerance*tolerance)
		return 1;
	return 0;
}
