#include <stdlib.h>
#include <stdint.h>

#include "nb/memory_bot.h"
#include "nb/interpolation_bot/nonpolynomial.h"

#define POW2(a) ((a)*(a))

static void eval_gi(uint8_t N, uint8_t dim, const double *ni,
		   const double *ri, const double *x,
		   double *eval, double (*g)(double));
static double eval_g(double x, double (*g)(double));
static double get_ri(const double *ri, uint8_t i);
static void eval_phi(uint8_t N, const double *gi, double *phi);
static double get_k(uint8_t N, const double *phi);

void nb_nonpolynomial_simple_eval(uint8_t N, uint8_t dim, const double *ni,
				  const double *x, double *eval)
{
	nb_nonpolynomial_custom_eval(N, dim, ni, NULL, x, eval, NULL);
}

void nb_nonpolynomial_simple_eval_deriv_x(uint8_t N, uint8_t dim,
					  const double *ni, const double *x,
					  double *eval);
void nb_nonpolynomial_simple_eval_deriv_y(uint8_t N, uint8_t dim,
					  const double *ni, const double *x,
					  double *eval);

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

	for (uint8_t i = 0; i < N; i++)
		eval[i] = phi[i] / k;

	NB_SOFT_FREE(memsize, memblock);
}

static void eval_gi(uint8_t N, uint8_t dim, const double *ni,
		    const double *ri, const double *x,
		    double *eval, double (*g)(double))
{
	for (uint8_t i = 0; i < N; i++) {
		double dist2 = 0;
		for (uint8_t d = 0; d < dim; d++) {
			double xi = ni[i*dim+d];
			dist2 += POW2(x[d] - xi);
		}
		double zi = dist2 / get_ri(ri, i);
		eval[i] = eval_g(zi, g);
	}
}

static double eval_g(double x, double (*g)(double))
{
	double gx;
	if (NULL != g)
		gx = g(x);
	else
		gx = x;
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

void nb_nonpolynomial_custom_eval_deriv_x(uint8_t N, uint8_t dim, const double *ni,
					  /* NULL for all ri = 1*/ const double *ri,
					  const double *x, double *eval,
					  /* NULL for g(x) = x */
					  double (*g)(double),
					  /* NULL for g'(x) = 1 */
					  double (*dg)(double));
void nb_nonpolynomial_custom_eval_deriv_y(uint8_t N, uint8_t dim, const double *ni,
					  /* NULL for all ri = 1*/ const double *ri,
					  const double *x, double *eval,
					  /* NULL for g(x) = x */
					  double (*g)(double),
					  /* NULL for g'(x) = 1 */
					  double (*dg)(double));
