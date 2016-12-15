#ifndef __NB_SOLVER_BOT_MATRIX_TRIANGULAR_H__
#define __NB_SOLVER_BOT_MATRIX_TRIANGULAR_H__

#include <stdint.h>

/* Triangular solvers for complete matrix */
void nb_matrix_forward_solve(const double *const L,
			      const double *const b,
			      double* _x, uint32_t N);
void nb_matrix_backward_solve(const double *const U,
			       const double *const b,
			       double* _x, uint32_t N);

#endif
