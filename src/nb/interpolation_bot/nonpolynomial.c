#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include "nb/memory_bot.h"
#include "nb/interpolation_bot/nonpolynomial.h"

#define POW2(a) ((a)*(a))
#define MAX(a,b) (((a)>(b))?(a):(b))

static void eval_phis(uint32_t N, uint8_t dim, const double *ni,
		      const double *ri, const double *x,
		      double *eval);
static double get_dist(uint8_t dim, const double *x1, const double *x2);
static double get_ri(const double *ri, uint32_t i);
static double func_z(uint8_t dim, const double *xi,
		     double ri, const double *x);
static void func_grad_z(uint8_t dim, const double *xi,
			double ri, const double *x,
			double *grad);
static double func_phi(double x);
static double func_dphi(double x);
static void eval_Pis(uint32_t N, const double *phis, double *eval);
static double get_Pi(uint32_t N, const double *phis, uint32_t i);
static double get_k(uint32_t N, const double *Pis);
static void eval_interpolators(uint32_t N, const double *Pis,
			       double k, double *eval);
static void eval_grad_Pis(uint32_t N, uint8_t dim, const double *ni,
			 const double *ri,  const double *phis,
			 const double *x, double *eval);
static void eval_grad_Pi(uint32_t N, uint8_t dim, const double *ni,
			   const double *ri, const double *phi, 
			   const double *x, uint32_t i, double *eval);
static double get_Pij(uint32_t N, const double *phi, uint32_t i, uint32_t j);
static void get_grad_k(uint32_t N, uint8_t dim, const double *grad_gi,
		       double *grad_k);
static void eval_grad_interpolators(uint32_t N, uint8_t dim,
				    const double *gi, const double *grad_gi,
				    double k, const double *grad_k,
				    double *eval);

void nb_nonpolynomial_eval(uint32_t N, uint8_t dim, const double *ni,
			   /* NULL for all ri = 1*/ const double *ri,
			   const double *x, double *eval)
{
	uint32_t memsize = 2 * N * sizeof(double);
	char *memblock = NB_SOFT_MALLOC(memsize);
	
	double *phis = (void*) memblock;
	double *Pis = (void*) (memblock + N * sizeof(double));
	
	eval_phis(N, dim, ni, ri, x, phis);

	eval_Pis(N, phis, Pis);

	double k = get_k(N, Pis);

	eval_interpolators(N, Pis, k, eval);

	NB_SOFT_FREE(memsize, memblock);
}

static void eval_phis(uint32_t N, uint8_t dim, const double *ni,
		      const double *ri, const double *x,
		      double *eval)
{
	for (uint32_t i = 0; i < N; i++) {
		double r = get_ri(ri, i);
		double z = func_z(dim, &(ni[i*dim]), r, x);
		eval[i] = func_phi(z);
	}
}

static double func_z(uint8_t dim, const double *xi,
		     double ri, const double *x)
{
	double d = get_dist(dim, x, xi);
	return d / ri;
}

static void func_grad_z(uint8_t dim, const double *xi,
			double ri, const double *x,
			double *grad)
{
	double dist = get_dist(dim, x, xi);
	
	dist = MAX(dist, 1e-6);

	for (uint8_t d = 0; d < dim; d++)
		grad[d] = (x[d] - xi[d]) / (ri * dist);
}

static double get_dist(uint8_t dim, const double *x1, const double *x2)
{
	double sum = 0;
	for (uint8_t d = 0; d < dim; d++)
		sum += POW2(x1[d] - x2[d]);
	return sqrt(sum);
}

static double get_ri(const double *ri, uint32_t i)
{
	double r;
	if (NULL != ri)
		r = ri[i];
	else
		r = 1.0;
	return r;
}

static double func_phi(double x)
{
	return x * log(x+1);
}

static double func_deriv_phi(double x)
{
	return log(x+1) + x/(x+1);
}

static void eval_Pis(uint32_t N, const double *phis, double *eval)
{
	for (uint32_t i = 0; i < N; i++)
		eval[i] = get_Pi(N, phis, i);
}

static double get_Pi(uint32_t N, const double *phis, uint32_t i)
{
	double prod = 1.0;
	for (uint32_t k = 0; k < N; k++) {
		if (k != i)
			prod *= phis[k];
	}
	return prod;
}

static double get_k(uint32_t N, const double *Pis)
{
	double k = 0;
	for (uint32_t i = 0; i < N; i++)
		k += Pis[i];
	return k;
}

static void eval_interpolators(uint32_t N, const double *Pis,
			       double k, double *eval)
{
	for (uint32_t i = 0; i < N; i++)
		eval[i] = Pis[i] / k;
}

void nb_nonpolynomial_eval_grad(uint32_t N, uint8_t dim,
				const double *ni,
				/* NULL for all ri = 1*/
				const double *ri, const double *x,
				double *eval)
{
	uint32_t memsize = ((2 + dim) * N + dim) * sizeof(double);
	char *memblock = NB_SOFT_MALLOC(memsize);
	
	double *phis = (void*) memblock;
	double *Pis = (void*) (memblock + N * sizeof(double));
	double *grad_Pis = (void*) (memblock + 2 * N * sizeof(double));
	double *grad_k = (void*)(memblock + (2 + dim) * N * sizeof(double));

	eval_phis(N, dim, ni, ri, x, phis);

	eval_Pis(N, phis, Pis);

	double k = get_k(N, Pis);

	eval_grad_Pis(N, dim, ni, ri, phis, x, grad_Pis);

	get_grad_k(N, dim, grad_Pis, grad_k);

	eval_grad_interpolators(N, dim, Pis, grad_Pis, k, grad_k, eval);

	NB_SOFT_FREE(memsize, memblock);	
}

static void eval_grad_Pis(uint32_t N, uint8_t dim, const double *ni,
			 const double *ri,  const double *phis,
			 const double *x, double *eval)
{
	for (uint32_t i = 0; i < N; i++)
		eval_grad_Pi(N, dim, ni, ri, phis, x, i, &(eval[i * dim]));
}

static void eval_grad_Pi(uint32_t N, uint8_t dim, const double *ni,
			 const double *ri, const double *phis, 
			 const double *x, uint32_t i, double *eval)
{
	uint32_t memsize = dim * sizeof(double);
	double *grad_z = NB_SOFT_MALLOC(memsize);
	memset(eval, 0, dim * sizeof(*eval));
	for (uint32_t j = 0; j < N; j++) {
		if (j != i) {			  
			double Pij = get_Pij(N, phis, i, j);
			double rj = get_ri(ri, j);
			double zj = func_z(dim, &(ni[j*dim]), rj, x);
			double dphi = func_deriv_phi(zj);
			func_grad_z(dim, &(ni[j*dim]), rj, x, grad_z);
			double prod = Pij * dphi;
			for (uint8_t d = 0; d < dim; d++) {
				eval[d] += prod * grad_z[d];
			}
		}
	}
	NB_SOFT_FREE(memsize, grad_z);
}

static double get_Pij(uint32_t N, const double *phis, uint32_t i, uint32_t j)
{
	double prod = 1.0;
	for (uint32_t k = 0; k < N; k++) {
		if (k != i && k != j)
			prod *= phis[k];
	}
	return prod;
}

static void get_grad_k(uint32_t N, uint8_t dim, const double *grad_Pis,
		       double *grad_k)
{
	memset(grad_k, 0, dim * sizeof(*grad_k));
	for (uint32_t i = 0; i < N; i++) {
		for (uint8_t d = 0; d < dim; d++)
			grad_k[d] += grad_Pis[i * dim + d];
	}
}

static void eval_grad_interpolators(uint32_t N, uint8_t dim,
				    const double *Pis, const double *grad_Pis,
				    double k, const double *grad_k,
				    double *eval)
{
	double k2 = POW2(k);
	for (uint32_t i = 0; i < N; i++) {
		for (uint8_t d = 0; d < dim; d++) {
			double dPi = grad_Pis[i * dim + d];
			double dk = grad_k[d];
			double divisor = (dPi * k - Pis[i] * dk);
			eval[i * dim + d] = divisor / k2;
		}
	}
}
