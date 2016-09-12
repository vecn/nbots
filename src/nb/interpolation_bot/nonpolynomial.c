#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include "nb/memory_bot.h"
#include "nb/interpolation_bot/nonpolynomial.h"

#define POW2(a) ((a)*(a))

static double thin_plate(double x);
static double d_thin_plate(double x);
static void eval_gi(uint8_t N, uint8_t dim, const double *ni,
		   const double *ri, const double *x,
		   double *eval, double (*g)(double));
static double get_dist(uint8_t dim, const double *x1, const double *x2);
static double eval_g(double x, double (*g)(double));
static double get_ri(const double *ri, uint8_t i);
static void eval_phi(uint8_t N, const double *gi, double *phi);
static double get_k(uint8_t N, const double *phi);
static void eval_interpolator(uint8_t N, const double *phi,
			      double k, double *eval);
static void eval_gi_and_dgi(uint8_t N, uint8_t dim, const double *ni,
			    const double *ri, const double *x,
			    double *gi, double *dgi,
			    double (*g)(double), double (*dg)(double));
static double eval_dg(double x, double (*dg)(double));
static void eval_grad_phi(uint8_t N, uint8_t dim, const double *ni,
			  const double *ri, const double *gi,
			  const double *dgi, const double *x,
			  double *grad_phi);
static double get_prod(uint8_t N, const double *gi, uint8_t i, uint8_t j);
static void get_grad_k(uint8_t N, uint8_t dim, const double *grad_phi,
		       double *grad_k);
static void eval_grad_interpolator(uint8_t N, uint8_t dim,
				   const double *phi, const double *grad_phi,
				   double k, const double *grad_k,
				   double *eval);

void nb_nonpolynomial_simple_eval(uint8_t N, uint8_t dim, const double *ni,
				  const double *x, double *eval)
{
	nb_nonpolynomial_custom_eval(N, dim, ni, NULL, x, eval, NULL);
}

void nb_nonpolynomial_simple_eval_grad(uint8_t N, uint8_t dim,
				       const double *ni, const double *x,
				       double *eval)
{
	nb_nonpolynomial_custom_eval_grad(N, dim, ni, NULL, x,
					  eval, NULL, NULL);
}

void nb_nonpolynomial_thin_plate_eval(uint8_t N, uint8_t dim, const double *ni,
				      const double *x, double *eval)
{
	nb_nonpolynomial_custom_eval(N, dim, ni, NULL, x, eval, thin_plate);
}

void nb_nonpolynomial_thin_plate_eval_grad(uint8_t N, uint8_t dim,
					   const double *ni, const double *x,
					   double *eval)
{
	nb_nonpolynomial_custom_eval_grad(N, dim, ni, NULL, x, eval,
					  thin_plate, d_thin_plate);
}

static double thin_plate(double x)
{
	return POW2(x) * log(x);
}

static double d_thin_plate(double x)
{
	return x * (1 + 2 * log(x));
}

void nb_nonpolynomial_custom_eval(uint8_t N, uint8_t dim, const double *ni,
				  /* NULL for all ri = 1*/ const double *ri,
				  const double *x, double *eval,
				  /* NULL for g(x) = x */ double (*g)(double))
{
	uint16_t memsize = 2 * N * sizeof(double);
	char *memblock = NB_SOFT_MALLOC(memsize);
	
	double *phi = (void*) memblock;
	double *gi = (void*) (memblock + N * sizeof(double));
	
	eval_gi(N, dim, ni, ri, x, gi, g);

	eval_phi(N, gi, phi);

	double k = get_k(N, phi);

	eval_interpolator(N, phi, k, eval);

	NB_SOFT_FREE(memsize, memblock);
}

static void eval_gi(uint8_t N, uint8_t dim, const double *ni,
		    const double *ri, const double *x,
		    double *eval, double (*g)(double))
{
	for (uint8_t i = 0; i < N; i++) {
		double di = get_dist(dim, x, &(ni[i*dim]));
		double zi = di / get_ri(ri, i);
		eval[i] = eval_g(zi, g);
	}
}

static double get_dist(uint8_t dim, const double *x1, const double *x2)
{
	double sum = 0;
	for (uint8_t i = 0; i < dim; i++)
		sum += POW2(x1[i] - x2[i]);
	return sqrt(sum);
}

static double eval_g(double x, double (*g)(double))
{
	double gx;
	if (NULL != g)
		gx = g(x);
	else
		gx = 1-exp(-x) + sqrt(x);/* AQUI voy */
	return gx;
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

static void eval_phi(uint8_t N, const double *gi, double *phi)
{
	for (uint8_t i = 0; i < N; i++) {
		phi[i] = 1.0;
		for (uint8_t j = 0; j < N; j++) {
			if (i != j)
				phi[i] *= gi[j];
		}
	}
}

static double get_k(uint8_t N, const double *phi)
{
	double k = 0;
	for (uint8_t i = 0; i < N; i++)
		k += phi[i];
	return k;
}

static void eval_interpolator(uint8_t N, const double *phi,
			      double k, double *eval)
{
	for (uint8_t i = 0; i < N; i++)
		eval[i] = phi[i] / k;
}

void nb_nonpolynomial_custom_eval_grad(uint8_t N, uint8_t dim,
				       const double *ni,
				       /* NULL for all ri = 1*/
				       const double *ri,
				       const double *x, double *eval,
				       /* NULL for g(x) = x */
				       double (*g)(double),
				       /* NULL for g'(x) = 1 */
				       double (*dg)(double))
{
	uint32_t memsize = ((3 + dim) * N + dim) * sizeof(double);
	char *memblock = NB_SOFT_MALLOC(memsize);
	
	double *phi = (void*) memblock;
	double *gi = (void*) (memblock + N * sizeof(double));
	double *dgi = (void*) (memblock + 2 * N * sizeof(double));
	double *grad_phi = (void*) (memblock + 3 * N * sizeof(double));
	double *grad_k = (void*)(memblock + (3 + dim) * N * sizeof(double));
	
	eval_gi_and_dgi(N, dim, ni, ri, x, gi, dgi, g, dg);

	eval_phi(N, gi, phi);

	eval_grad_phi(N, dim, ni, ri, gi, dgi, x, grad_phi);

	double k = get_k(N, phi);

	get_grad_k(N, dim, grad_phi, grad_k);

	eval_grad_interpolator(N, dim, phi, grad_phi, k, grad_k, eval);

	NB_SOFT_FREE(memsize, memblock);	
}

static void eval_gi_and_dgi(uint8_t N, uint8_t dim, const double *ni,
			    const double *ri, const double *x,
			    double *gi, double *dgi,
			    double (*g)(double), double (*dg)(double))
{
	for (uint8_t i = 0; i < N; i++) {
		double di = get_dist(dim, x, &(ni[i*dim]));
		double zi = di / get_ri(ri, i);
		gi[i] = eval_g(zi, g);
		dgi[i] = eval_dg(zi, dg);
	}
}

static double eval_dg(double x, double (*dg)(double))
{
	double dgx;
	if (NULL != dg)
		dgx = dg(x);
	else
		dgx = exp(-x) + 1/(2*sqrt(x));
	return dgx;
}

static void eval_grad_phi(uint8_t N, uint8_t dim, const double *ni,
			  const double *ri, const double *gi,
			  const double *dgi, const double *x,
			  double *grad_phi)
{
	memset(grad_phi, 0, dim * N * sizeof(*grad_phi));
	for (uint8_t i = 0; i < N; i++) {
		for (uint8_t j = 0; j < N; j++) {
			if (j != i) {			  
				double c1 = get_prod(N, gi, i, j);
				double rj = get_ri(ri, j);
				double dj = get_dist(dim, x, &(ni[j*dim]));
				double c2 = dgi[j] / (dj * rj);
				for (uint8_t d = 0; d < dim; d++) {
					double h = x[d] - ni[j * dim + d];
					double eval = c1 * c2 * h;
					grad_phi[i * dim + d] += eval;
				}
			}
		}
	}
}

static double get_prod(uint8_t N, const double *gi, uint8_t i, uint8_t j)
{
	double prod = 1.0;
	for (uint8_t k = 0; k < N; k++) {
		if (k != i && k != j)
			prod *= gi[k];
	}
	return prod;
}

static void get_grad_k(uint8_t N, uint8_t dim, const double *grad_phi,
		       double *grad_k)
{
	memset(grad_k, 0, dim * sizeof(*grad_k));
	for (uint8_t i = 0; i < N; i++) {
		for (uint8_t d = 0; d < dim; d++)
			grad_k[d] += grad_phi[i * dim + d];
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
