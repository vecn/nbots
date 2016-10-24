#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "nb/memory_bot.h"
#include "nb/solver_bot/sparse/sparse.h"
#include "nb/solver_bot/sparse/eigen/lanczos.h"

#include "../sparse_struct.h"

#define POW2(a) ((a)*(a))

void nb_sparse_eigen_lanczos(const nb_sparse_t* const A,
			      double *_eigenmax,/* Out */ 
			      double *_eigenmin,/* Out */
			      int* it,          /* Out */
			      double tolerance,
			      uint32_t omp_parallel_threads)
{
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
	int i, j;

	/* Declare structures and variables to be used */
	double* alpha = nb_allocate_zero_mem(A->N * sizeof(double));
	double* beta = nb_allocate_zero_mem(A->N * sizeof(double));
	double* v = nb_allocate_zero_mem(A->N * sizeof(double));
	double* w = nb_allocate_zero_mem(A->N * sizeof(double));
	double* v_prev = nb_allocate_zero_mem(A->N * sizeof(double));

	*it = 0;
	beta[*it] = 0;
	v[0] = 1;
	/* Start iterations */
	double normw2 = 1;
	double delta_eigen = 1;

	while (delta_eigen >= POW2(tolerance) &&
	       normw2 >= tolerance && *it < A->N) {
		/* Step 1 and 2 */
		double ak = 0;
#pragma omp parallel for reduction(+:ak) num_threads(omp_parallel_threads) private(i, j)
		for (i = 0; i < A->N; i++) {
			w[i] = 0;
			for (j = 0; j < A->rows_size[i]; j++)
				w[i] += A->rows_values[i][j]*v[A->rows_index[i][j]];
			w[i] -= beta[*it]*v_prev[i];
			ak += w[i]*v[i];
		}
		alpha[*it] = ak;
		/* Step 3 and 4 */
		*it = *it+1;
		normw2 = 0;
		for (i = 0; i < A->N; i++) {
			w[i] -= ak*v[i];
			normw2 += w[i]*w[i];
		}
		normw2 = sqrt(normw2);
		beta[*it] = normw2;
		/* Step 5 */
		for (i = 0; i < A->N; i++) {
			v_prev[i] = v[i];
			v[i] = w[i]/beta[*it];
		}
		/* Step 6 and 7 */
		if (*it > 1) {
			double delta_eig1 = *_eigenmax;
			double delta_eigk = *_eigenmin;
			nb_sparse_eigen_givens(alpha, beta, 1, _eigenmax, tolerance, *it);
			nb_sparse_eigen_givens(alpha, beta, *it, _eigenmin, tolerance, *it);
      
			delta_eigen = fabs(delta_eig1-*_eigenmax);
			double tmp = fabs(delta_eigk-*_eigenmin);
			if (tmp > delta_eigen)
				delta_eigen = tmp;
		} else {
			*_eigenmax = ak;
			*_eigenmin = ak;
		}
	}

	/* Free memory */
	nb_free_mem(alpha);
	nb_free_mem(beta);
	nb_free_mem(w);
	nb_free_mem(v);
	nb_free_mem(v_prev);
}
