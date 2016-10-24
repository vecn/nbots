#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "nb/memory_bot.h"
#include "nb/solver_bot/sparse/sparse.h"
#include "nb/solver_bot/sparse/eigen/givens.h"

#include "../sparse_struct.h"

void nb_sparse_eigen_givens(const double* const main_diag, 
			     const double* const uplw_diag,
			     int i, double *_eigenvalue,
			     double tolerance, uint32_t N)
{
	/* Compute the eigen values of a symmetric tridiagonal 
	 * matrix (using Givens method), where
	 *   > *main_diag is the main diagonal of the matrix.
	 *   > *uplw_diag is the upper and lower diagonal, the
	 *     first value will be ignored.
	 *   > 'i' is the position of the eigen value to compute, 
	 *     assuming:
	 *        l1 > l2 > l2 > ... > ln
	 *     where li is the ith eigenvalue.
	 *   > The program return through '*_eigenvalue' the
	 *     resulting value.
	 */
	/* Iterative variables */
	int l;
	/* Intitialize and allocate structures */
	double* p = nb_allocate_zero_mem((N+1) * sizeof(double));

	/* Init algorithm */
	int k = N;
	double a = main_diag[0]-fabs(uplw_diag[1]);
	double tmp = main_diag[k-1]-fabs(uplw_diag[k-1]);
	if (tmp < a)
		a = tmp;
	double b = main_diag[0]+fabs(uplw_diag[1]);
	tmp = main_diag[k-1]+fabs(uplw_diag[k-1]);
	if (tmp > b)
		b = tmp;
	for (l=1; l<k-1; l++) {
		tmp = main_diag[l]-
			fabs(uplw_diag[l+1])-fabs(uplw_diag[l]);
		if (tmp < a)
			a = tmp;
		tmp = main_diag[l]+
			fabs(uplw_diag[l+1])+fabs(uplw_diag[l]);
		if (tmp > b)
			b = tmp;
	}
	/* Init iterations */
	while (fabs(b-a) > (fabs(a)+fabs(b)) * tolerance) {
		/* Step 1 */
		*_eigenvalue = (a+b)/2;
		/* Step 2 */
		double r = 0;
		double s = 0;
		p[0] = 1;
		p[1] = main_diag[0] - *_eigenvalue;
		for (l = 1; l < k; l++) {
			p[l+1] = (main_diag[l]-*_eigenvalue)*p[l]-
				(uplw_diag[l])*(uplw_diag[l])*p[l-1];
		}
		for (l=1; l<k+1; l++) {
			if (p[l]*p[l-1] <= 0)
				r++;
			if (p[l] == 0)
				s++;
		}
		double gamma = r-s;
		/* Step 3 */
		if (gamma > k-i)
			b = *_eigenvalue;
		else
			a = *_eigenvalue;
	}

	/* Free memory */
	nb_free_mem(p);
}
