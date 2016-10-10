#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "nb/math_bot.h"
#include "nb/memory_bot.h"
#include "nb/container_bot/array.h"
#include "nb/solver_bot/matrix.h"

#include "sparse_struct.h"

#define POW2(a) ((a)*(a))

static inline void matrix_r_solver(const double *const A, int n, 
				   double *d, double *b);

double vcn_matrix_2X2_inverse(const double *const A, double* A_inv)
{
	double det = (A[0] * A[3] - A[1] * A[2]);
	A_inv[0] = A[3] / det;
	A_inv[1] = -A[1] / det;
	A_inv[2] = -A[2] / det;
	A_inv[3] = A[0] / det;
	return det;
}

void vcn_matrix_2X2_eigen(const double *const A,
			  double* Lambda, /* Output (diag) */
			  double* P,      /* Output */
			  double tolerance)
/* Returns A decomposed into P Lambda P' */
{
	if (fabs(A[1]) < tolerance && fabs(A[2]) < tolerance) {
		Lambda[0] = A[0];
		Lambda[1] = A[3];
		P[0] = 1.0;
		P[2] = 0.0;
		P[1] = 0.0;
		P[3] = 1.0;
	} else {
		double T = A[0] + A[3];
		double D = A[0] * A[3] - A[1] * A[2];
		double root = sqrt(POW2(T)/4.0 - D);
		Lambda[0] = T/2.0 + root;
		Lambda[1] = T/2.0 - root;
		if (fabs(A[2]) > tolerance) {
			P[0] = Lambda[0] - A[3];
			P[2] = A[2];
			P[1] = Lambda[1] - A[3];
			P[3] = A[2];
		} else {
			P[0] = A[1];
			P[2] = Lambda[0] - A[0];
			P[1] = A[1];
			P[3] = Lambda[1] - A[0];
		}
	}
	if (fabs(Lambda[0]) < fabs(Lambda[1])) {
		double aux = Lambda[0];
		Lambda[0] = Lambda[1];
		Lambda[1] = aux;

		aux = P[0];
		P[0] = P[1];
		P[1] = aux;

		aux = P[2];
		P[2] = P[3];
		P[3] = aux;
	}
	/* Normalize eigenvectors */
	double normalizer = sqrt(POW2(P[0]) + POW2(P[2]));
	P[0] /= normalizer;
	P[2] /= normalizer;
	normalizer = sqrt(POW2(P[1]) + POW2(P[3]));
	P[1] /= normalizer;
	P[3] /= normalizer;
}
double vcn_matrix_2X2_det(double *A)
{
	return A[0] * A[3] - A[1] * A[2];
}

double vcn_matrix_3X3_det(double *A)
{
	return 
		A[0] * A[4] * A[8] + 
		A[3] * A[7] * A[2] +
		A[6] * A[1] * A[5] -
		A[2] * A[4] * A[6] -
		A[5] * A[7] * A[0] -
		A[8] * A[1] * A[3];
}

double vcn_matrix_2X2_inverse_destructive(double *A)
{
	double det = A[0] * A[3] - A[1] * A[2];
	double tmp = A[0];
	A[0] = A[3]/det;
	A[3] = tmp/det;
	A[1] = -A[1]/det;
	A[2] = -A[2]/det;
	return det;
}

double vcn_matrix_3X3_inverse_destructive(double *A)
{
	double det = vcn_matrix_3X3_det(A);
	double a11, a12, a13, a21, a22, a23, a31, a32, a33;
	a11 = A[0];
	a12 = A[1];
	a13 = A[2];
	a21 = A[3];
	a22 = A[4];
	a23 = A[5];
	a31 = A[6];
	a32 = A[7];
	a33 = A[8];
	A[0] = (a33*a22-a32*a23)/det;
	A[1] = -(a33*a12-a32*a13)/det;
	A[2] = (a23*a12-a22*a13)/det;
	A[3] = -(a33*a21-a31*a23)/det;
	A[4] = (a33*a11-a31*a13)/det;
	A[5] = -(a23*a11-a21*a13)/det;
	A[6] = (a32*a21-a31*a22)/det;
	A[7] = -(a32*a11*a31*a12)/det;
	A[8] = (a22*a11-a21*a12)/det;
	return det;
}

int vcn_matrix_cholesky_decomposition(const double *const A,
				      double* _LplusLt,   /* Out */
				      uint32_t N)
{
	/* The 'L' computed correspond to [L + Lt] */
	for (uint32_t j=0; j < N; j++) {
		_LplusLt[j*N+j] = A[j*N+j];
		for (uint32_t k=0; k<j; k++)
			_LplusLt[j*N+j] -= _LplusLt[j*N+k]*_LplusLt[j*N+k];
		if (_LplusLt[j*N+j] <= 1e-16)
			return 1;
		_LplusLt[j*N+j] = sqrt(_LplusLt[j*N+j]);
		for (uint32_t i=j+1;  i<N; i++) {
			_LplusLt[i*N+j] = A[i*N+j];
			for (uint32_t k=0; k<j; k++)
				_LplusLt[j*N+j] -= _LplusLt[i*N+k]*_LplusLt[j*N+k];
			_LplusLt[i*N+j] /= _LplusLt[j*N+j];
			_LplusLt[j*N+i] = _LplusLt[i*N+j];
		}
	}
	return 0;
}

void vcn_matrix_cholesky_solve
(const double *const LplusLt,
 const double *const b, 
 double* _x,            /* Out */
 uint32_t N)
/* Solve the system LL'x = b, where LL'= A */
{
	double* z = nb_allocate_zero_mem(N * sizeof(double));
	vcn_matrix_forward_solve(LplusLt, b, z, N);
	vcn_matrix_backward_solve(LplusLt, z, _x, N);
	/* Free memory */
	nb_free_mem(z);
}

double vcn_matrix_cond1(const double *const A, int n){
	double *Acopy = (double*)nb_allocate_mem(n*n*sizeof(double));
	memcpy(Acopy, A, n*n*sizeof(double));
	double *x = (double*)nb_allocate_mem(n*sizeof(double));
	double *p = (double*)nb_allocate_mem(n*sizeof(double));
	double *pm = (double*)nb_allocate_mem(n*sizeof(double));
	/* Compute QR decomposition */
	int sing;
	double *c = (double*)nb_allocate_mem(n*sizeof(double));
	double *diag = (double*)nb_allocate_mem(n*sizeof(double));
	vcn_matrix_qr_decomposition(Acopy, n, c, diag, &sing);
	/* Compute ||A||_1 (Max absolute column sum) */
	double estimation = fabs(diag[0]);
	for(uint32_t i=1; i<n; i++){
		double sum = 0;
		for(uint32_t j=0; j<i; j++)
			sum +=  fabs(Acopy[i*n+j]);
		double tmp = fabs(diag[i]) + sum;
		estimation = (tmp>estimation)?tmp:estimation;
	}
	/* Solve R^T x = e, selecting e as they proceed */
	x[0] = 1.0/diag[0];
	for(uint32_t i=1; i<n; i++)
		p[i] = Acopy[i]*x[0];
	for(uint32_t i=1; i<n; i++){
		/* Select ej and calculate xj */
		double xp = (1-p[i])/diag[i];
		double xm = (-1-p[i])/diag[i];
		double tmp = fabs(xp);
		double tmpm = fabs(xm);
		for(uint32_t j=i+1; j<n; j++){
			pm[j] = p[j] + Acopy[i*n+j]*xm;
			tmpm += (fabs(pm[j])/fabs(diag[j]));
			p[j] = p[j] + Acopy[i*n+j]*xp;
			tmp += (fabs(p[j])/fabs(diag[j]));
		}
		if(tmp > tmpm) x[i] = xp;    /* ej = 1 */
		else{
			/* ej = -1 */
			x[i] = xm;
			for(uint32_t j=i+1; j<n; j++)
				p[j]  = pm[j];
		}
	}
	double xnorm = 0;
	for(uint32_t i=0; i<n; i++) xnorm += fabs(x[i]);
	estimation /= xnorm;
	/* Solve Ry = x */
	vcn_matrix_qr_solve(Acopy, n, c, diag, x);

	xnorm = 0;
	for(uint32_t i=0; i<n; i++) xnorm += fabs(x[i]);
	estimation *= xnorm;
	/* Free memory */
	nb_free_mem(Acopy);
	nb_free_mem(x);
	nb_free_mem(p);
	nb_free_mem(pm);
	nb_free_mem(diag);
	nb_free_mem(c);
	return estimation;
}

double vcn_matrix_cond2(const double *const A, int n){
	double *Acopy = (double*)nb_allocate_mem(n*n*sizeof(double));
	memcpy(Acopy, A, n*n*sizeof(double));
	double* w = (double*)nb_allocate_mem(n*sizeof(double));
	double* V = (double*)nb_allocate_mem(n*n*sizeof(double));
	/* Compute SVD Decomposition */ 
	vcn_matrix_svd_decomposition(Acopy, w, V, n, n);
	/* Compute condition number estimation */
	double estimation = w[0]/w[n-1];
	/* Free memory */
	nb_free_mem(Acopy);
	nb_free_mem(w);
	nb_free_mem(V);
	return estimation;
}

void vcn_matrix_qr_decomposition(double *A, /* Overwritten */
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

void vcn_matrix_qr_solve(const double *const A,
			 int n, double *c, double *d,
			 double *b /* Solution overwritten */){
	/* Numerical Recipes in C
	 * Solves the set of n linear equations Ax = b. A, c and d are the 
	 * input as the output of the routine vcn_matrix_qr_decomposition  and
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

void vcn_matrix_svd_decomposition(double *A, /* Overwritten wit U */
				  double *w, double *V, 
				  int n, int m){
	/* Numerical Recipes in C.
	 * Given  a  matrix  A,  this  routine  computes its singular value 
	 * decomposition, A = UWV^T. The matrix U replaces A on output. The
	 * diagonal  matrix  of  singular  values W is output as a vector w. 
	 * The matrix V (not the transpose V^T) is output as V.
	 */
	int flag, i, its, j, jj, k, l, nm;
	double anorm, c, f, g, h, s, scale, x, y, z, *rv1;
	rv1 = (double*)nb_allocate_zero_mem(n * sizeof(double));
	g = scale = anorm = 0.0;
	for(i=0; i<n; i++){
		l = i+1;
		rv1[i] = scale*g;
		g = s = scale = 0.0;
		if(i < m){
			for(k=i; k<m; k++) scale += fabs(A[k*n+i]);
			if(scale){
				for(k=i; k<m; k++){
					A[k*n+i] /= scale;
					s += A[k*n+i]*A[k*n+i];
				}
				f = A[i*n+i];
				g = -(f >= 0.0 ? fabs(sqrt(s)) : -fabs(sqrt(s)));
				h = f*g-s;
				A[i*n+i] = f-g;
				for(j=l; j<n; j++){
					for(s=0.0, k=i; k<m; k++) s += A[k*n+i]*A[k*n+j];
					f = s/h;
					for(k=i; k<m; k++) A[k*n+j] += f*A[k*n+i];
				}
				for(k=i; k<m; k++) A[k*n+i] *= scale;
			}
		}
		w[i] = scale*g;
		g = s = scale = 0.0;
		if(i < m && i != n-1){
			for(k=l; k<n; k++) scale += fabs(A[i*n+k]);
			if(scale){
				for(k=l; k<n; k++){
					A[i*n+k] /= scale;
					s += A[i*n+k]*A[i*n+k];
				}
				f = A[i*n+l];
				g = -(f >= 0.0 ? fabs(sqrt(s)) : -fabs(sqrt(s)));
				h = f*g-s;
				A[i*n+l] = f-g;
				for(k=l; k<n; k++) rv1[k] = A[i*n+k]/h;
				for(j=l; j<m; j++){
					for(s = 0.0, k=l; k<n; k++) s += A[j*n+k]*A[i*n+k];
					for(k=l; k<n; k++) A[j*n+k] += s*rv1[k];
				}
				for(k=l; k<n; k++) A[i*n+k] *= scale;
			}
		}
		anorm = 
			(anorm >= (fabs(w[i])+fabs(rv1[i]))) ? anorm: (fabs(w[i])+fabs(rv1[i]));
	}
	for(i=n-1; i>=0; i--){
		if(i<n-1){
			if(g){
				for(j=l; j<n; j++)
					V[j*n+i] = (A[i*n+j]/A[i*n+l])/g;
				for(j=l; j<n; j++){
					for(s=0.0, k=l; k<n; k++) s += A[i*n+k]*V[k*n+j];
					for(k=l; k<n; k++) V[k*n+j] += s*V[k*n+i];
				}
			}
			for(j=l; j<n; j++) V[i*n+j] = V[j*n+i] = 0.0;
		}
		V[i*n+i] = 1.0;
		g = rv1[i];
		l = i;
	}
	for(i=((n<m)?(n-1):(m-1)); i>=0; i--){
		l = i+1;
		g = w[i];
		for(j=l; j<n; j++) A[i*n+j] = 0.0;
		if(g){
			g = 1.0/g;
			for(j=l; j<n; j++){
				for(s=0.0, k=l; k<m; k++) s += A[k*n+i]*A[k*n+j];
				f = (s/A[i*n+i])*g;
				for(k=i; k<m; k++) A[k*n+j] += f*A[k*n+i];
			}
			for(j=i; j<m; j++) A[j*n+i] *= g;
      
		}else for(j=i; j<n; j++) A[j*n+i] = 0.0;
		++A[i*n+i];
	}
	for(k=n-1; k>=0; k--){
		for(its=0; its<30; its++){
			flag = 1;
			for(l=k; l>=0; l--){
				nm = l-1;
				if(fabs(rv1[l])+anorm == anorm){
					flag = 0;
					break;
				}
				if(fabs(w[nm])+anorm == anorm) break;
			}
			if(flag){
				c = 0.0;
				s = 1.0;
				for(i=l; i<k; i++){
					f = s*rv1[i];
					rv1[i] = c*rv1[i];
					if(fabs(f)+anorm == anorm) break;
					g = w[i];
					h = vcn_math_hypo(f,g);
					w[i] = h;
					h = 1.0/h;
					c = g*h;
					s = -f*h;
					for(j=0; j<m; j++){
						y = A[j*n+nm];
						z = A[j*n+i];
						A[j*n+nm] = y*c+z*s;
						A[j*n+i] = z*c-y*s;
					}
				}
			}
			z = w[k];
			if(l == k){
				if(z < 0.0){
					w[k] = -z;
					for(j=0; j<n; j++) 
						V[j*n+k] = -V[j*n+k];
				}
				break;
			}
			if(its == 29){
				printf("SVD Error 1\n");
				exit(1);
			}
			x = w[l];
			nm = k-1;
			y = w[nm];
			g = rv1[nm];
			h = rv1[k];
			f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g = vcn_math_hypo(f, 1.0);
			f = ((x-z)*(x+z) + 
			     h*((y/(f + (f >= 0.0 ? fabs(g):-fabs(g))))-h))/x;
			c = s = 1.0;
			for(j=l; j<=nm; j++){
				i = j+1;
				g = rv1[i];
				y = w[i];
				h = s*g;
				g = c*g;
				z = vcn_math_hypo(f, h);
				rv1[j] = z;
				c = f/z;
				s = h/z;
				f = x*c+g*s;
				g = g*c-x*s;
				h = y*s;
				y *= c;
				for(jj=0; jj<n; jj++){
					x = V[jj*n+j];
					z = V[jj*n+i];
					V[jj*n+j] = x*c+z*s;
					V[jj*n+i] = z*c-x*s;
				}
				z = vcn_math_hypo(f, h);
				w[j] = z;
				if(z){
					z = 1.0/z;
					c = f*z;
					s = h*z;
				}
				f = c*g+s*y;
				x = c*y-s*g;
				for(jj=0; jj<m; jj++){
					y = A[jj*n+j];
					z = A[jj*n+i];
					A[jj*n+j] = y*c+z*s;
					A[jj*n+i] = z*c-y*s;
				}
			}
			rv1[l] = 0.0;
			rv1[k] = f;
			w[k] = x;
		}
	}
	/* Free memory */
	nb_free_mem(rv1);
}


void vcn_matrix_svd_solve(const double *const U,
			  const double *const w, 
			  const double *const V, 
			  double *x,             /* Out */
			  const double *const b, 
			  int n, int m){
	/* Numerical Recipes in C.
	 * Solves Ax = b for a vector X, where A is specified by the arrays u, w, v 
	 * as returned by svd_dcmp, m and n are the dimensions of a, and will be equal
	 * for square matrices. b is the input right-hand side. x is the output solution 
	 * vector. No input quantities are destroyed, so the routine may be called 
	 * sequentially with different b's. 
	 */
	int jj, j, i;
	double s, *tmp;
	tmp = (double*)nb_allocate_mem(n*sizeof(double));
	for(j=0; j<n; j++){
		s = 0.0;
		if(w[j]){
			for(i=0; i<m; i++) s += U[i*n+j]*b[i];
			s /= w[j];
		}
		tmp[j] = s;
	}
	for(j=0; j<n; j++){
		s = 0.0;
		for(jj=0; jj<n; jj++) s += V[j*n+jj]*tmp[jj];
		x[j] = s;
	}
	nb_free_mem(tmp);
}


void vcn_matrix_forward_solve(const double *const L,
			      const double *const b, 
			      double *_x, uint32_t N)
/* Solve the system Lx = b, where L is a lower triangular matrix */
{
	/* This solver cannot be parallelized to guaranty the solution */
	for(uint32_t i=0; i< N; i++){
		_x[i] = b[i];
		for(uint32_t k=0; k < i; k++)
			_x[i] -= L[i*N+k]*_x[k];
		_x[i] /= L[i*N+i];
	}
}

void vcn_matrix_backward_solve(const double *const U,
			       const double *const b,
			       double *_x, uint32_t N)
/* Solve the system Ux = b, where U is a upper triangular matrix */
{
	/* This solver cannot be parallelized to guaranty the solution */
	for(int i=N-1; i>=0 ; i--){
		_x[i] = b[i];
		for(int k=i+1; k<N; k++)
			_x[i] -= U[i*N+k]*_x[k];
		_x[i] /= U[i*N+i];
	}
}
