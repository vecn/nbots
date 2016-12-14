#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "nb/memory_bot.h"
#include "nb/math_bot.h"
#include "nb/solver_bot/matrix/svd.h"

void nb_matrix_svd_decomposition(double *A, /* Overwritten wit U */
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
					h = nb_math_hypo(f,g);
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
			g = nb_math_hypo(f, 1.0);
			f = ((x-z)*(x+z) + 
			     h*((y/(f + (f >= 0.0 ? fabs(g):-fabs(g))))-h))/x;
			c = s = 1.0;
			for(j=l; j<=nm; j++){
				i = j+1;
				g = rv1[i];
				y = w[i];
				h = s*g;
				g = c*g;
				z = nb_math_hypo(f, h);
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
				z = nb_math_hypo(f, h);
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


void nb_matrix_svd_solve(const double *const U,
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
