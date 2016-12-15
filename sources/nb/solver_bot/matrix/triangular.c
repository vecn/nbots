#include <stdint.h>

#include "nb/solver_bot/matrix/triangular.h"

void nb_matrix_forward_solve(const double *const L,
			      const double *const b, 
			      double *_x, uint32_t N)
/* Solve the system Lx = b, where L is a lower triangular matrix */
{
	for (uint32_t i = 0; i < N; i++) {
		_x[i] = b[i];
		for (uint32_t k = 0; k < i; k++)
			_x[i] -= L[i*N+k]*_x[k];
		_x[i] /= L[i*N+i];
	}
}

void nb_matrix_backward_solve(const double *const U,
			       const double *const b,
			       double *_x, uint32_t N)
/* Solve the system Ux = b, where U is a upper triangular matrix */
{
	for (int32_t i = N-1; i >= 0 ; i--) {
		_x[i] = b[i];
		for (int32_t k = i+1; k < N; k++)
			_x[i] -= U[i*N+k]*_x[k];
		_x[i] /= U[i*N+i];
	}
}
