#ifndef __NB_SOLVER_BOT_MATRIX_TRIDIAGONAL_SOLVER_H__
#define __NB_SOLVER_BOT_MATRIX_TRIDIAGONAL_SOLVER_H__

void nb_matrix_solve_tridiagonal_destructive(uint32_t N, const double *ld,
					     const double *md, double *ud,
					     const double *b, double *x);
void nb_matrix_solve_tridiagonal(uint32_t N, const double *ld,
				 const double *md, const double *ud,
				 const double *b, double *x);

#endif
