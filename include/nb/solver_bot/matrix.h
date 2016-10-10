#ifndef __NB_SOLVER_BOT_MATRIX_H__
#define __NB_SOLVER_BOT_MATRIX_H__

#include <stdint.h>

double vcn_matrix_2X2_inverse(const double *const A, double* A_inv);
  
void vcn_matrix_2X2_eigen(const double *const A,
			  /* Return: A decomposed into P Lambda P' */
			  double* Lambda,  /* Output (diag) */
			  double* P,       /* Output */
			  double tolerance);
double vcn_matrix_2X2_det(double *A);

double vcn_matrix_3X3_det(double *A);

double vcn_matrix_2X2_inverse_destructive(double *A);

double vcn_matrix_3X3_inverse_destructive(double *A);

int vcn_matrix_cholesky_decomposition(const double *const A,
				      double* _LplusLt,      /* Out */
				      uint32_t N);
void vcn_matrix_cholesky_solve(const double *const LplusLt,
			       const double *const b, 
			       double* _x,            /* Out */
			       uint32_t N);
double vcn_matrix_cond1(const double *const A, int N);
double vcn_matrix_cond2(const double *const A, int N);
void vcn_matrix_qr_decomposition(double * A, /* Overwritten */
				 int N, 
				 double *c, double *d, int *sing);
void vcn_matrix_qr_solve(const double *const A,
			 int N, double *c, double *d,
			 double *b /* Solution overwritten */);
void vcn_matrix_svd_decomposition(double *A, /* Overwritten with U */
				  double *w, double *V, 
				  int N, int M);
void vcn_matrix_svd_solve(const double *const U,
			  const double *const w, 
			  const double *const V, 
			  double *x,             /* Out */
			  const double *const b, int N, int M);

/* Triangular solvers for complete matrix */
void vcn_matrix_forward_solve(const double *const L,
			      const double *const b,
			      double* _x, uint32_t N);
void vcn_matrix_backward_solve(const double *const U,
			       const double *const b,
			       double* _x, uint32_t N);

#endif
