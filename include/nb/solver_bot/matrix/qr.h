#ifndef __NB_SOLVER_BOT_MATRIX_QR_H__
#define __NB_SOLVER_BOT_MATRIX_QR_H__

void nb_matrix_qr_decomposition(double * A, /* Overwritten */
				 int N, 
				 double *c, double *d, int *sing);
void nb_matrix_qr_solve(const double *const A,
			 int N, double *c, double *d,
			 double *b /* Solution overwritten */);

#endif
