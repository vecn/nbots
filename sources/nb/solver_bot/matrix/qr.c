#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "nb/memory_bot.h"
#include "nb/solver_bot/matrix/qr.h"

static inline void matrix_r_solver(const double *const A, int n, 
				   double *d, double *b);

void nb_matrix_qr_decomposition(double *A, /* Overwritten */
				 int n,
				 double *c, double *d, int *sing){
	/* Numerical Recipes in C.
	 * Constructs the QR decomposition of A. The upper triangular matrix
	 * R is returned in the upper triangle of A, except for the diagonal 
	 * elements of R which are returned in d[1..n]. The orthogonal matrix 
	 * Q is represented as a product of n-1 Householder matrices Q_1...Q_{n-1},
	 * where Q_j = 1 - u_j prod u_j/c_j. The ith component of u_j is zero
	 * for i=1,...,j-1 while the nonzero components are returned in A_ij 
	 * for i=j,...,n.  "sing" return as true(1) if singularity is encountered
	 * during the decomposition, but the decomposition is still completed in 
	 * this case, otherwise it returns false(0).
	 */
	double scale, sigma, sum, tau;
	*sing = 0;
	for(int k=0; k<n-1; k++){
		scale = 0.0;
		for(int i=k; i<n; i++) 
			scale = (scale > fabs(A[i*n+k]))?(scale):(fabs(A[i*n+k]));
		if(scale == 0.0){
			*sing = 1;
			c[k] = d[k] = 0.0;
		}else{
			for(int i=k; i<n; i++)
				A[i*n+k] /= scale;
			sum = 0.0;
			for(int i=k; i<n; i++) 
				sum += A[i*n+k]*A[i*n+k];
			sigma = (A[k*n+k] >= 0.0 ? fabs(sqrt(sum)) : -fabs(sqrt(sum)));
			A[k*n+k] += sigma;
			c[k] = sigma*A[k*n+k];
			d[k] = -scale*sigma;
			for(int j=k+1; j<n; j++){
				sum = 0.0;
				for(int i=k; i<n; i++)
					sum += A[i*n+k]*A[i*n+j];
				tau = sum/c[k];
				for(int i=k; i<n; i++)
					A[i*n+j] -= tau*A[i*n+k];
			}
		}
	}
	d[n-1] = A[n*n-1];
	if(d[n-1] == 0.0) *sing = 1;
}

static inline void matrix_r_solver(const double *const A, int n, double *d, double *b){
	b[n-1] /= d[n-1];
	for(int i=n-2; i>=0; i--){
		double sum = 0.0;
		for(int j=i+1; j<n; j++) sum += A[i*n+j]*b[j];
		b[i] = (b[i]-sum)/d[i];
	}
}

void nb_matrix_qr_solve(const double *const A,
			 int n, double *c, double *d,
			 double *b /* Solution overwritten */){
	/* Numerical Recipes in C
	 * Solves the set of n linear equations Ax = b. A, c and d are the 
	 * input as the output of the routine nb_matrix_qr_decomposition  and
	 * are not modified. b is input as the right hand side vector, and 
	 * is overwritten with the solution vector on output.
	 */
	for(int j=0; j<n-1; j++){
		double sum = 0.0;
		for(int i=j; i<n; i++)
			sum += A[i*n+j]*b[i];
		double tau = sum/c[j];
		for(int i=j; i<n; i++)
			b[i] -= tau*A[i*n+j];
	}
	matrix_r_solver(A, n, d, b);
}
