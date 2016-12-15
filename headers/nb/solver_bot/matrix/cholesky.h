#ifndef __NB_SOLVER_BOT_MATRIX_CHOLESKY_H__
#define __NB_SOLVER_BOT_MATRIX_CHOLESKY_H__

int nb_matrix_cholesky_decomposition(const double *const A,
				      double* _LplusLt,      /* Out */
				      uint32_t N);
void nb_matrix_cholesky_solve(const double *const LplusLt,
			       const double *const b, 
			       double* _x,            /* Out */
			       uint32_t N);

#endif
