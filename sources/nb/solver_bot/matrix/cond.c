#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "nb/memory_bot.h"
#include "nb/solver_bot/matrix/qr.h"
#include "nb/solver_bot/matrix/svd.h"
#include "nb/solver_bot/matrix/cond.h"

double nb_matrix_cond1(const double *const A, int n){
	double *Acopy = (double*)nb_allocate_mem(n*n*sizeof(double));
	memcpy(Acopy, A, n*n*sizeof(double));
	double *x = (double*)nb_allocate_mem(n*sizeof(double));
	double *p = (double*)nb_allocate_mem(n*sizeof(double));
	double *pm = (double*)nb_allocate_mem(n*sizeof(double));
	/* Compute QR decomposition */
	int sing;
	double *c = (double*)nb_allocate_mem(n*sizeof(double));
	double *diag = (double*)nb_allocate_mem(n*sizeof(double));
	nb_matrix_qr_decomposition(Acopy, n, c, diag, &sing);
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
	nb_matrix_qr_solve(Acopy, n, c, diag, x);

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

double nb_matrix_cond2(const double *const A, int n){
	double *Acopy = (double*)nb_allocate_mem(n*n*sizeof(double));
	memcpy(Acopy, A, n*n*sizeof(double));
	double* w = (double*)nb_allocate_mem(n*sizeof(double));
	double* V = (double*)nb_allocate_mem(n*n*sizeof(double));
	/* Compute SVD Decomposition */ 
	nb_matrix_svd_decomposition(Acopy, w, V, n, n);
	/* Compute condition number estimation */
	double estimation = w[0]/w[n-1];
	/* Free memory */
	nb_free_mem(Acopy);
	nb_free_mem(w);
	nb_free_mem(V);
	return estimation;
}
