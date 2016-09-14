#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include "nb/memory_bot.h"
#include "nb/interpolation_bot/nonpolynomial.h"

#define POW2(a) ((a)*(a))

static void eval_phi(uint8_t N, uint8_t dim, const double *ni,
		   const double *ri, const double *x,
		   double *eval);
static double get_dist2(uint8_t dim, const double *x1, const double *x2);
static double get_ri(const double *ri, uint8_t i);
static void eval_gi(uint8_t N, const double *phi, double c, double *eval);
static double get_Pi(uint8_t N, const double *phi, uint8_t i);
static double get_k(uint8_t N, const double *gi);
static void eval_interpolator(uint8_t N, const double *gi,
			      double k, double *eval);
static void eval_grad_gi(uint8_t N, uint8_t dim, const double *ni,
			 const double *ri, const double *phi,
			 const double *x, double c, double *eval);
static double eval_grad_Pi(uint8_t N, uint8_t dim, const double *ni,
			   const double *ri, const double *phi, 
			   const double *x, uint8_t i, double *eval);
static double get_Pij(uint8_t N, const double *phi, uint8_t i, uint8_t j);
static void get_grad_k(uint8_t N, uint8_t dim, const double *grad_gi,
		       double *grad_k);
static void eval_grad_interpolator(uint8_t N, uint8_t dim,
				   const double *gi, const double *grad_gi,
				   double k, const double *grad_k,
				   double *eval);

void nb_nonpolynomial_eval(uint8_t N, uint8_t dim, const double *ni,
			   /* NULL for all ri = 1*/ const double *ri,
			   const double *x, double c, double *eval)
{
	uint32_t memsize = 2 * N * sizeof(double);
	char *memblock = NB_SOFT_MALLOC(memsize);
	
	double *phi = (void*) memblock;
	double *gi = (void*) (memblock + N * sizeof(double));
	
	eval_phi(N, dim, ni, ri, x, phi);

	eval_gi(N, phi, c, gi);

	double k = get_k(N, gi);

	eval_interpolator(N, gi, k, eval);

	NB_SOFT_FREE(memsize, memblock);
}

static void eval_phi(uint8_t N, uint8_t dim, const double *ni,
		    const double *ri, const double *x,
		    double *eval)
{
	for (uint8_t i = 0; i < N; i++) {
		double di2 = get_dist2(dim, x, &(ni[i*dim]));
		double r = get_ri(ri, i);
		double zi2 = di2 / POW2(r);
		eval[i] = zi2;
	}
}

static double get_dist2(uint8_t dim, const double *x1, const double *x2)
{
	double sum = 0;
	for (uint8_t i = 0; i < dim; i++)
		sum += POW2(x1[i] - x2[i]);
	return sum;
}

static double get_ri(const double *ri, uint8_t i)
{
	double r;
	if (NULL != ri)
		r = ri[i];
	else
		r = 1.0;
	return r;
}

static void eval_gi(uint8_t N, const double *phi, double c, double *eval)
{
	for (uint8_t i = 0; i < N; i++)
		eval[i] = 2 * get_Pi(N, phi, i) / (1 + pow(phi[i], c));
	for (uint8_t i = 0; i < N; i++)/* TEMPORAL */
		eval[i] = get_Pi(N, phi, i);
}

static double get_Pi(uint8_t N, const double *phi, uint8_t i)
{
	double prod = 1.0;
	for (uint8_t k = 0; k < N; k++) {
		if (k != i)
			prod *= phi[k];
	}
	return prod;
}

static double get_k(uint8_t N, const double *gi)
{
	double k = 0;
	for (uint8_t i = 0; i < N; i++)
		k += gi[i];
	return k;
}

static void eval_interpolator(uint8_t N, const double *gi,
			      double k, double *eval)
{
	for (uint8_t i = 0; i < N; i++)
		eval[i] = gi[i] / k;
}

void nb_nonpolynomial_eval_grad(uint8_t N, uint8_t dim,
				const double *ni,
				/* NULL for all ri = 1*/
				const double *ri, const double *x,
				double c, double *eval)
{
	uint32_t memsize = ((2 + dim) * N + dim) * sizeof(double);
	char *memblock = NB_SOFT_MALLOC(memsize);
	
	double *phi = (void*) memblock;
	double *gi = (void*) (memblock + N * sizeof(double));
	double *grad_gi = (void*) (memblock + 2 * N * sizeof(double));
	double *grad_k = (void*)(memblock + (2 + dim) * N * sizeof(double));

	eval_phi(N, dim, ni, ri, x, phi);

	eval_gi(N, phi, c, gi);

	eval_grad_gi(N, dim, ni, ri, phi, x, c, grad_gi);

	double k = get_k(N, gi);

	get_grad_k(N, dim, grad_gi, grad_k);

	eval_grad_interpolator(N, dim, gi, grad_gi, k, grad_k, eval);

	NB_SOFT_FREE(memsize, memblock);	
}

static void eval_grad_gi(uint8_t N, uint8_t dim, const double *ni,
			 const double *ri,  const double *phi,
			 const double *x, double c, double *eval)
{
	double *grad_Pi = alloca(dim * sizeof(*grad_Pi));
	for (uint8_t i = 0; i < N; i++) {
		double Pi = eval_grad_Pi(N, dim, ni, ri, phi, x, i, grad_Pi);
		double phib = pow(phi[i], c - 1);
		double phic = pow(phi[i], c);
		double div = POW2(1 + phic);
		double subs = Pi * c * phib;
		for (uint8_t d = 0; d < dim; d++) {
			double sum = grad_Pi[d] * (1 + phic);
			double prod = sum - subs;
			eval[i*dim + d] = 2 * prod / div;
			eval[i*dim + d] = grad_Pi[d];/* TEMPORAL */
		}
	}
}

static double eval_grad_Pi(uint8_t N, uint8_t dim, const double *ni,
			 const double *ri, const double *phi, 
			   const double *x, uint8_t i, double *eval)
{
	double Pi = 1.0;
	memset(eval, 0, dim * sizeof(*eval));
	for (uint8_t j = 0; j < N; j++) {
		if (j != i) {			  
			double rj = get_ri(ri, j);
			double Pij = get_Pij(N, phi, i, j);
			double prod = 2 * Pij / POW2(rj);
			for (uint8_t d = 0; d < dim; d++) {
				double h = x[d] - ni[j * dim + d];
				eval[d] += prod * h;
			}
			Pi = phi[j] * Pij; /* Same cheap operation */
		}
	}
	return Pi;
}

static double get_Pij(uint8_t N, const double *phi, uint8_t i, uint8_t j)
{
	double prod = 1.0;
	for (uint8_t k = 0; k < N; k++) {
		if (k != i && k != j)
			prod *= phi[k];
	}
	return prod;
}

static void get_grad_k(uint8_t N, uint8_t dim, const double *grad_gi,
		       double *grad_k)
{
	memset(grad_k, 0, dim * sizeof(*grad_k));
	for (uint8_t i = 0; i < N; i++) {
		for (uint8_t d = 0; d < dim; d++)
			grad_k[d] += grad_gi[i * dim + d];
	}
}

static void eval_grad_interpolator(uint8_t N, uint8_t dim,
				   const double *phi, const double *grad_phi,
				   double k, const double *grad_k,
				   double *eval)
{
	double k2 = POW2(k);
	for (uint8_t i = 0; i < N; i++) {
		for (uint8_t d = 0; d < dim; d++) {
			double dphi = grad_phi[i * dim + d];
			double dk = grad_k[d];
			double divisor = (dphi * k - phi[i] * dk);
			eval[i * dim + d] = divisor / k2;
		}
	}
}
