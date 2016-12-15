#ifndef __NB_SOLVER_BOT_MATRIX_SVD_H__
#define __NB_SOLVER_BOT_MATRIX_SVD_H__

void nb_matrix_svd_decomposition(double *A, /* Overwritten with U */
				  double *w, double *V, 
				  int N, int M);
void nb_matrix_svd_solve(const double *const U,
			  const double *const w, 
			  const double *const V, 
			  double *x,             /* Out */
			  const double *const b, int N, int M);


#endif
