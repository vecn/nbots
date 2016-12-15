#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "nb/memory_bot.h"
#include "nb/solver_bot/matrix/triangular.h"
#include "nb/solver_bot/matrix/cholesky.h"

int nb_matrix_cholesky_decomposition(const double *const A,
				      double* _LplusLt,   /* Out */
				      uint32_t N)
{
	/* The 'L' computed correspond to [L + Lt] */
	for (uint32_t j=0; j < N; j++) {
		_LplusLt[j*N+j] = A[j*N+j];
		for (uint32_t k=0; k<j; k++)
			_LplusLt[j*N+j] -= _LplusLt[j*N+k]*_LplusLt[j*N+k];
		if (_LplusLt[j*N+j] <= 1e-16)
			return 1;
		_LplusLt[j*N+j] = sqrt(_LplusLt[j*N+j]);
		for (uint32_t i=j+1;  i<N; i++) {
			_LplusLt[i*N+j] = A[i*N+j];
			for (uint32_t k=0; k<j; k++)
				_LplusLt[j*N+j] -= _LplusLt[i*N+k]*_LplusLt[j*N+k];
			_LplusLt[i*N+j] /= _LplusLt[j*N+j];
			_LplusLt[j*N+i] = _LplusLt[i*N+j];
		}
	}
	return 0;
}

void nb_matrix_cholesky_solve
(const double *const LplusLt,
 const double *const b, 
 double* _x,            /* Out */
 uint32_t N)
/* Solve the system LL'x = b, where LL'= A */
{
	double* z = nb_allocate_zero_mem(N * sizeof(double));
	nb_matrix_forward_solve(LplusLt, b, z, N);
	nb_matrix_backward_solve(LplusLt, z, _x, N);
	/* Free memory */
	nb_free_mem(z);
}
