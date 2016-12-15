#ifndef __NB_SOLVER_BOT_MATRIX_MATRIX2X2_H__
#define __NB_SOLVER_BOT_MATRIX_MATRIX2X2_H__

double nb_matrix_2X2_inverse(const double *const A, double* A_inv);
  
void nb_matrix_2X2_eigen(const double *const A,
			  /* Return: A decomposed into P Lambda P' */
			  double* Lambda,  /* Output (diag) */
			  double* P,       /* Output */
			  double tolerance);
double nb_matrix_2X2_det(double *A);
double nb_matrix_2X2_inverse_destructive(double *A);

#endif
