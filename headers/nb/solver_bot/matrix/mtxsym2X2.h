#ifndef __NB_SOLVER_BOT_MATRIX_MTXSYM2X2_H__
#define __NB_SOLVER_BOT_MATRIX_MTXSYM2X2_H__

double nb_mtxsym_2X2_inverse(const double A[3], double A_inv[3]);

void nb_mtxsym_2X2_eigen(const double A[3], double Lambda[2],
			 double P[4], double tol);
double nb_mtxsym_2X2_det(const double A[3]);
double nb_mtxsym_2X2_inverse_destructive(double A[3]);

#endif
