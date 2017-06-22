#ifndef __NB_SOLVER_BOT_MATRIX_MATRIX2X2_H__
#define __NB_SOLVER_BOT_MATRIX_MATRIX2X2_H__

double nb_matrix_2X2_inverse(const double A[4], double A_inv[4]);
  
void nb_matrix_2X2_eigen(const double A[4],
			  /* Return: A decomposed into P Lambda P' */
			  double* Lambda,  /* Output (diag) */
			  double* P,       /* Output */
			  double tolerance);
double nb_matrix_2X2_det(const double A[4]);
double nb_matrix_2X2_inverse_destructive(double A[4]);

#endif
