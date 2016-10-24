#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "nb/memory_bot.h"
#include "nb/solver_bot/sparse/sparse.h"
#include "nb/solver_bot/sparse/eigen/inv_power.h"

#include "../sparse_struct.h"

#define POW2(a) ((a)*(a))

int nb_sparse_eigen_ipower(const nb_sparse_t *const A,
			    nb_solver_t solver,
			    int h, double mu, 
			    double **_eigenvecs,/* Out */ 
			    double *_eigenvals, /* Out */
			    int* it,            /* Out */
			    double tolerance,
			    uint32_t omp_parallel_threads)
{
	/* The program must receive all the pointers allocated, where
	 *  > A is a nb_sparse_t matrix
	 *  > _eigenvecs is an array of size h to store h eigenvectors.
	 *  > _eigenvals is an array of size h to store the h
	 *    eigenvalues approximated.
	 *  > h is the number of eigenvalues to be computed.
	 *  > '*it' will store (after computation) the iterations needed 
	 *    to compute each eigenvalue (is a return value).
	 */

	/* Declare structures and variables to be used */
	uint32_t i, j, c, d; /* Iterative variables */
	double pnorm, rnorm2;
	/* Allocate memory for structures */
	double* p = nb_allocate_zero_mem(A->N * sizeof(double));
	double* z = nb_allocate_zero_mem(A->N * sizeof(double));

	/* Set M in A (copy ptr to modify const A, it will be restored) */
	nb_sparse_t* A_ptr_copy = (nb_sparse_t*)A;
	if (mu != 0.0)
		for (c = 0; c < A->N; c++)
			nb_sparse_add(A_ptr_copy, c, c, -mu);            /* M = A - mu*I */

	/* LU Decomposition in case of LU Solver */
	nb_sparse_t *L = NULL;
	nb_sparse_t *U = NULL;
	if (NB_SOLVER_CHK == solver) {
		nb_sparse_alloc_LU(A, &L, &U);
		nb_sparse_decompose_Cholesky(A, L, U, omp_parallel_threads);
	} else if (NB_SOLVER_LUD == solver) {
		nb_sparse_alloc_LU(A, &L, &U);
		nb_sparse_decompose_LU(A, L, U, omp_parallel_threads);
	}

	/* Deflation inverse power method */
	for (i = 0; i < h; i++) {
		it[i] = 0;
		rnorm2 = 1;
		/* Initialize q0 such that ||qk||=1 */
		_eigenvecs[i][0] = 1;
		for(j=1; j< A->N; j++)
			_eigenvecs[i][j] = 0;

		/* Start loop */
		double rnorm2_diff = 1;
		while (rnorm2 > POW2(tolerance) &&
		       rnorm2_diff > POW2(tolerance)) {
			/* Step 1 */
			if (NB_SOLVER_CHK == solver ||
			    NB_SOLVER_LUD == solver)
				nb_sparse_solve_LU(L, U, _eigenvecs[i], p);
			else if (NB_SOLVER_CGJ == solver)
				nb_sparse_solve_CG_precond_Jacobi(A,_eigenvecs[i], p, 
								   nb_sparse_get_size(A)*10, 1e-3, 
								   NULL, NULL,
								   omp_parallel_threads);
			else
				nb_sparse_solve_conjugate_gradient(A,_eigenvecs[i], p, 
								    nb_sparse_get_size(A)*10, 1e-3, 
								    NULL, NULL,
								    omp_parallel_threads);
			/* Step 2 */
			pnorm = nb_vector_get_norm(p, A->N);
			for(c=0; c < A->N; c++)
				_eigenvecs[i][c] = p[c]/pnorm;
			/* Step 3 */
			for (j = 0; j < i; j++){
				double alpha = 0;

#pragma omp parallel for reduction(+:alpha) num_threads(omp_parallel_threads) schedule(guided) private(c)
				for(c=0; c < A->N; c++)
					alpha += _eigenvecs[i][c]*_eigenvecs[j][c];

#pragma omp parallel for num_threads(omp_parallel_threads) private(c)
				for(c=0; c < A->N; c++)
					_eigenvecs[i][c] -= alpha*_eigenvecs[j][c];
			}
			/* Step 4 */
			/* Paralelize the operation zk = A*qk */
#pragma omp parallel for num_threads(omp_parallel_threads) schedule(guided) private(d) private(c)
			for(c=0; c < A->N; c++){
				z[c] = 0;
				for(d=0; d < A->rows_size[c]; d++){
					double aii = A->rows_values[c][d];
					if(c == d) aii += mu;
					z[c] += aii
						*(_eigenvecs[i][A->rows_index[c][d]]);
				}
			}
			/* Step 5 */
			double sigma = 0;
			for(c=0; c < A->N; c++)
				sigma += _eigenvecs[i][c]*z[c];
			_eigenvals[i] = sigma;
			/* Step 6 and 7 */
			rnorm2_diff = rnorm2;
			rnorm2 = 0;
			for(c=0; c < A->N; c++)
				rnorm2 += POW2(z[c]-sigma*_eigenvecs[i][c]);
			rnorm2_diff = fabs(rnorm2_diff-rnorm2);
			it[i]++;
		}
	}
	/* Restore A */
	if (mu != 0.0)
		for (c = 0; c < A->N; c++)
			nb_sparse_add(A_ptr_copy, c, c, mu);  /* A = M + mu*I */
    
	/* Destroy LU decomposition */
	if (NB_SOLVER_CHK == solver || NB_SOLVER_LUD == solver) {
		nb_sparse_destroy(U);
		nb_sparse_destroy(L);
	}
	/* Free memory */
	nb_free_mem(p);
	nb_free_mem(z);

	return 0;
}
