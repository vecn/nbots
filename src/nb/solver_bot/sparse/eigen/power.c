#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "nb/memory_bot.h"
#include "nb/solver_bot/sparse/sparse.h"
#include "nb/solver_bot/sparse/solvers/gauss_seidel.h"

#include "../sparse_struct.h"

#define POW2(a) ((a)*(a))

void nb_sparse_eigen_power(const nb_sparse_t* const A, int h,
			    double **_eigenvecs,/* Out */ 
			    double *_eigenvals, /* Out */
			    int *it,            /* Out */
			    double tolerance,
			    uint32_t omp_parallel_threads){
	/* The program must receive all the pointers allocated, where
	 *  > A is a nb_sparse_t matrix
	 *  > _eigenvecs is an array of size h to store h eigenvectors.
	 *  > _eigenvals is an array of size h to store the h greatest
	 *    eigenvalues approximated.
	 *  > h is the number of eigenvalues to be computed.
	 *  > '*it' will store (after computation) the iterations needed 
	 *    to compute each eigenvalue (is a return value).
	 */
		
	/* Declare structures and variables to be used */
	uint32_t i, j, c, d; /* Iterative variables */
	double pnorm, rnorm2;
	/* Allocate memory for structures */
	double *p = nb_allocate_zero_mem(A->N * sizeof(double));

	/* Deflation power method */
	for (i = 0; i < h; i++) {
		it[i] = 0;
		rnorm2 = 1;
		/* Initialize q0 such that ||qk||=1 */
		_eigenvecs[i][0] = 1;
		for (j = 1; j < A->N; j++)
			_eigenvecs[i][j] = 0;

		for (c = 0; c < A->N; c++) {
			p[c] = 0;
			if(A->rows_index[c][0] == 0)
				p[c] = A->rows_values[c][0];
		}
		/* Start loop */
		while (rnorm2 > POW2(tolerance)) {
			/* Step 1 */
			pnorm = nb_vector_get_norm(p, A->N);
			for (c = 0; c < A->N; c++)
				_eigenvecs[i][c] = p[c]/pnorm;
			/* Step 2 */
			for (j = 0; j < i; j++) {
				double alpha = 0;

#pragma omp parallel for reduction(+:alpha) num_threads(omp_parallel_threads) schedule(guided) private(c)
				for(c=0; c < A->N; c++)
					alpha += _eigenvecs[i][c]*_eigenvecs[j][c];

#pragma omp parallel for num_threads(omp_parallel_threads) private(c)
				for(c=0; c < A->N; c++)
					_eigenvecs[i][c] -= alpha*_eigenvecs[j][c];
			}
			/* Step 3 */
			/* Paralelize the operation pk = A*qk */

#pragma omp parallel for schedule(guided) num_threads(omp_parallel_threads) private(c, d)
			for(c=0; c < A->N; c++){
				p[c] = 0;
				for(d=0; d < A->rows_size[c]; d++)
					p[c] += A->rows_values[c][d]
						*_eigenvecs[i][A->rows_index[c][d]];
			}
			/* Step 4 */
			double lambda = 0;
			for(c=0; c < A->N; c++)
				lambda += _eigenvecs[i][c]*p[c];
			_eigenvals[i] = lambda;
			/* Step 5 and 6 */
			rnorm2 = 0;
			for(c=0; c < A->N; c++)
				rnorm2 += POW2(p[c]-lambda*_eigenvecs[i][c]);
			it[i]++;
		}
	}

	/* Free memory */
	nb_free_mem(p);
}
