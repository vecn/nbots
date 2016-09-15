#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include "nb/memory_bot.h"
#include "nb/interpolation_bot/nonpolynomial.h"

#define POW2(a) ((a)*(a))
#define MAX(a,b) (((a)>(b))?(a):(b))

static void eval_phis(uint8_t N, uint8_t dim, const double *ni,
		      const double *ri, const double *x,
		      double *eval);
static double get_dist(uint8_t dim, const double *x1, const double *x2);
static double get_ri(const double *ri, uint8_t i);
static double func_z(uint8_t dim, const double *xi,
		     double ri, const double *x);
static void func_grad_z(uint8_t dim, const double *xi,
			double ri, const double *x,
			double *grad);
static double func_phi(double x);
static double func_dphi(double x);
static void eval_gs(uint8_t N, const double *phis, double c, double *eval);
static double get_Pi(uint8_t N, const double *phis, uint8_t i);
static double get_k(uint8_t N, const double *gs);
static void eval_interpolators(uint8_t N, const double *gs,
			       double k, double *eval);
static void eval_grad_gs(uint8_t N, uint8_t dim, const double *ni,
			 const double *ri,  const double *phis,
			 const double *x, double c, double *eval);
static double eval_grad_Pi(uint8_t N, uint8_t dim, const double *ni,
			   const double *ri, const double *phi, 
			   const double *x, uint8_t i, double *eval);
static double get_Pij(uint8_t N, const double *phi, uint8_t i, uint8_t j);
static void get_grad_k(uint8_t N, uint8_t dim, const double *grad_gi,
		       double *grad_k);
static void eval_grad_interpolators(uint8_t N, uint8_t dim,
				    const double *gi, const double *grad_gi,
				    double k, const double *grad_k,
				    double *eval);

void nb_nonpolynomial_eval(uint8_t N, uint8_t dim, const double *ni,
			   /* NULL for all ri = 1*/ const double *ri,
			   const double *x, double c, double *eval)
{
	uint32_t memsize = 2 * N * sizeof(double);
	char *memblock = NB_SOFT_MALLOC(memsize);
	
	double *phis = (void*) memblock;
	double *gs = (void*) (memblock + N * sizeof(double));
	
	eval_phis(N, dim, ni, ri, x, phis);

	eval_gs(N, phis, c, gs);

	double k = get_k(N, gs);

	eval_interpolators(N, gs, k, eval);

	NB_SOFT_FREE(memsize, memblock);
}

static void eval_phis(uint8_t N, uint8_t dim, const double *ni,
		      const double *ri, const double *x,
		      double *eval)
{
	for (uint8_t i = 0; i < N; i++) {
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
	for (uint8_t i = 0; i < dim; i++)
		sum += POW2(x1[i] - x2[i]);
	return sqrt(sum);
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

static double func_phi(double x)
{
	return sqrt(x);
}

static double func_deriv_phi(double x)
{
	return 0.5 / sqrt(x);
}

static void eval_gs(uint8_t N, const double *phis, double c, double *eval)
{
	for (uint8_t i = 0; i < N; i++) {
		double phic = pow(phis[i], c);
		eval[i] = get_Pi(N, phis, i) / (1 + phic);
	}
}

static double get_Pi(uint8_t N, const double *phis, uint8_t i)
{
	double prod = 1.0;
	for (uint8_t k = 0; k < N; k++) {
		if (k != i)
			prod *= phis[k];
	}
	return prod;
}

static double get_k(uint8_t N, const double *gs)
{
	double k = 0;
	for (uint8_t i = 0; i < N; i++)
		k += gs[i];
	return k;
}

static void eval_interpolators(uint8_t N, const double *gs,
			       double k, double *eval)
{
	for (uint8_t i = 0; i < N; i++)
		eval[i] = gs[i] / k;
}

void nb_nonpolynomial_eval_grad(uint8_t N, uint8_t dim,
				const double *ni,
				/* NULL for all ri = 1*/
				const double *ri, const double *x,
				double c, double *eval)
{
	uint32_t memsize = ((2 + dim) * N + dim) * sizeof(double);
	char *memblock = NB_SOFT_MALLOC(memsize);
	
	double *phis = (void*) memblock;
	double *gs = (void*) (memblock + N * sizeof(double));
	double *grad_gs = (void*) (memblock + 2 * N * sizeof(double));
	double *grad_k = (void*)(memblock + (2 + dim) * N * sizeof(double));

	eval_phis(N, dim, ni, ri, x, phis);

	eval_gs(N, phis, c, gs);

	eval_grad_gs(N, dim, ni, ri, phis, x, c, grad_gs);

	double k = get_k(N, gs);

	get_grad_k(N, dim, grad_gs, grad_k);

	eval_grad_interpolators(N, dim, gs, grad_gs, k, grad_k, eval);

	NB_SOFT_FREE(memsize, memblock);	
}

static void eval_grad_gs(uint8_t N, uint8_t dim, const double *ni,
			 const double *ri,  const double *phis,
			 const double *x, double c, double *eval)
{
	uint32_t memsize = 2 * dim * sizeof(double);
	char *memblock = NB_SOFT_MALLOC(memsize);
	double *grad_Pi = (void*) memblock;
	double *grad_z = (void*) (memblock + dim * sizeof(double));
	for (uint8_t i = 0; i < N; i++) {
		double Pi = eval_grad_Pi(N, dim, ni, ri, phis, x, i, grad_Pi);
		double phicm1 = pow(phis[i], c - 1);
		double phic = pow(phis[i], c);
		double div = POW2(1 + phic);
		double r = get_ri(ri, i);
		double z = func_z(dim, &(ni[i*dim]), r, x);
		double dphi = func_deriv_phi(z);
		double subs_prod = Pi * c * phicm1 * dphi;
		func_grad_z(dim, &(ni[i*dim]), r, x, grad_z);
		for (uint8_t d = 0; d < dim; d++) {
			double sum = grad_Pi[d] * (1 + phic);
			double subs = subs_prod * grad_z[d];
			double prod = sum - subs;
			eval[i*dim + d] = prod / div;
		}
	}
	NB_SOFT_FREE(memsize, memblock);
}

static double eval_grad_Pi(uint8_t N, uint8_t dim, const double *ni,
			 const double *ri, const double *phis, 
			   const double *x, uint8_t i, double *eval)
{
	uint32_t memsize = dim * sizeof(double);
	double *grad_z = NB_SOFT_MALLOC(memsize);
	double Pi = 1.0;
	memset(eval, 0, dim * sizeof(*eval));
	for (uint8_t j = 0; j < N; j++) {
		if (j != i) {			  
			double Pij = get_Pij(N, phis, i, j);
			double rj = get_ri(ri, j);
			double z = func_z(dim, &(ni[j*dim]), rj, x);
			double dphi = func_deriv_phi(z);
			func_grad_z(dim, &(ni[j*dim]), rj, x, grad_z);
			double prod = Pij * dphi;
			for (uint8_t d = 0; d < dim; d++) {
				eval[d] += prod * grad_z[d];
			}
			Pi *= phis[j];
		}
	}
	NB_SOFT_FREE(memsize, grad_z);
	return Pi;
}

static double get_Pij(uint8_t N, const double *phis, uint8_t i, uint8_t j)
{
	double prod = 1.0;
	for (uint8_t k = 0; k < N; k++) {
		if (k != i && k != j)
			prod *= phis[k];
	}
	return prod;
}

static void get_grad_k(uint8_t N, uint8_t dim, const double *grad_gs,
		       double *grad_k)
{
	memset(grad_k, 0, dim * sizeof(*grad_k));
	for (uint8_t i = 0; i < N; i++) {
		for (uint8_t d = 0; d < dim; d++)
			grad_k[d] += grad_gs[i * dim + d];
	}
}

static void eval_grad_interpolators(uint8_t N, uint8_t dim,
				    const double *gs, const double *grad_gs,
				    double k, const double *grad_k,
				    double *eval)
{
	double k2 = POW2(k);
	for (uint8_t i = 0; i < N; i++) {
		for (uint8_t d = 0; d < dim; d++) {
			double dgi = grad_gs[i * dim + d];
			double dk = grad_k[d];
			double divisor = (dgi * k - gs[i] * dk);
			eval[i * dim + d] = divisor / k2;
		}
	}
}
