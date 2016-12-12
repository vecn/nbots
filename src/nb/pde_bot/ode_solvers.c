#include <stdlib.h>
#include <stdint.h>

#include "nb/memory_bot.h"

#include "nb/pde_bot/ode_solvers.h"

#define POW2(a) ((a)*(a))

double eval_func(double (*f)(double x), double x);

void nb_fd_solve_ode_h2(uint32_t N, double L,
			double (*k)(double x), /* Diffusion (can be NULL) */
			double (*a)(double x), /* Convection (can be NULL) */
			double (*b)(double x), /* Absorption (can be NULL) */
			double (*c)(double x), /* Source (can be NULL) */
			double f0, double fL,  /* Boundary conditions */
			double *f_out)
/* -k f''(x) + a f'(x) + b f(x) = c */
/* h2: Error O(h2)                  */
{
	if (N < 3)
		goto EXIT;

	uint32_t n = N - 2;
	uint32_t memsize = 4 * n * sizeof(double);
	char *memblock = nb_soft_allocate_mem(memsize);
	double *ld = (void*) memblock;
	double *md = (void*) (memblock + n * sizeof(double));
	double *ud = (void*) (memblock + 2 * n * sizeof(double));
	double *rhs = (void*) (memblock + 3 * n * sizeof(double));

	double h = L/(N-1);
	for (uint32_t i = 0; i < n; i++) {
		double x = h + i * h;
		double kx = eval_func(k, x);
		double ax = eval_func(a, x);
		double bx = eval_func(b, x);
		double cx = eval_func(c, x);
		ld[i] = (-0.5 * ax * h - kx);
		md[i] = (bx * POW2(h) + 2.0 * kx);
		ud[i] = ( 0.5 * ax * h - kx);
		rhs[i] = POW2(h) * cx;
	}
	rhs[0] -= (-0.5 * eval_func(a, h) * h - eval_func(k, h)) * f0;
	rhs[n-1] -= (0.5 * eval_func(a, L - h) * h - eval_func(k, L - h)) * fL;

	f_out[0] = f0;
	f_out[N-1] = fL;
	nb_matrix_solve_tridiagonal_destructive(n, ld, md, ud, rhs,
						&(f_out[1]));

	nb_soft_free_mem(memsize, memblock);

EXIT:
	return;
}

double eval_func(double (*f)(double x), double x)
{
	double out;
	if (NULL != f)
		out = f(x);
	else
		out = 0.0;
	return out;
}

void nb_fem_solve_ode_o1(uint32_t N, double L,
			 double (*k)(double x), /* Diffusion (can be NULL) */
			 double (*a)(double x), /* Convection (can be NULL) */
			 double (*b)(double x), /* Absorption (can be NULL) */
			 double (*c)(double x), /* Source (can be NULL) */
			 double f0, double fL,  /* Boundary conditions */
			 double *f_out)
/* -k(x) f'' + a(x) f' + b(x) f = c(x) */
/* o1: Linear elements                 */
{
	if (N < 3)
		goto EXIT;

	uint32_t n = N - 2;
	uint32_t memsize = 4 * n * sizeof(double);
	char *memblock = nb_soft_allocate_mem(memsize);
	double *ld = (void*) memblock;
	double *md = (void*) (memblock + n * sizeof(double));
	double *ud = (void*) (memblock + 2 * n * sizeof(double));
	double *rhs = (void*) (memblock + 3 * n * sizeof(double));

	double h = L/(N-1);
	for (uint32_t i = 0; i < n; i++) {
		double x = h + i * h;
		double kx = eval_func(k, x);
		double ax = eval_func(a, x);
		double bx = eval_func(b, x);
		double cx = eval_func(c, x);
		ld[i] = (-0.5 * ax * h - kx + 0.16667 * bx * POW2(h));
		md[i] = (bx * POW2(h) + 2.0 * kx + 0.66666 * bx * POW2(h));
		ud[i] = ( 0.5 * ax * h - kx + 0.16667 * bx * POW2(h));
		rhs[i] = POW2(h) * cx;
	}
	rhs[0] -= (-0.5 * eval_func(a, h) * h - eval_func(k, h)
		   + 0.16667 * eval_func(b, h) * POW2(h)) * f0;
	rhs[n-1] -= (0.5 * eval_func(a, L - h) * h - eval_func(k, L - h)
		     + 0.16667 * eval_func(b, L - h) * POW2(h)) * fL;

	f_out[0] = f0;
	f_out[N-1] = fL;
	nb_matrix_solve_tridiagonal_destructive(n, ld, md, ud, rhs,
						&(f_out[1]));

	nb_soft_free_mem(memsize, memblock);

EXIT:
	return;
}

void nb_analytic_fd_solve_diffusion_convection
			(uint32_t N, double L,
			 double (*k)(double x), /* Diffusion (can be NULL) */
			 double (*a)(double x), /* Convection (can be NULL) */
			 double (*c)(double x), /* Source (can be NULL) */
			 double f0, double fL,  /* Boundary conditions */
			 double *f_out)
/* -k(x) f'' + a(x) f' = c(x) */
{
	if (N < 3)
		goto EXIT;

	uint32_t n = N - 2;
	uint32_t memsize = 4 * n * sizeof(double);
	char *memblock = nb_soft_allocate_mem(memsize);
	double *ld = (void*) memblock;
	double *md = (void*) (memblock + n * sizeof(double));
	double *ud = (void*) (memblock + 2 * n * sizeof(double));
	double *rhs = (void*) (memblock + 3 * n * sizeof(double));

	double h = L/(N-1);
	for (uint32_t i = 0; i < n; i++) {
		double x = h + i * h;
		double kx = eval_func(k, x);
		double ax = eval_func(a, x);
		double cx = eval_func(c, x);
		double p = ax / kx;
		double ep = 1.0 / (exp(p*h) + 1);
		ld[i] = ep - 1.0;
		md[i] = 1.0;
		ud[i] = -ep;
		rhs[i] = h * (cx/ax) * (1 - 2*ep);
	}
	double p0 = eval_func(a, h) / eval_func(k, h);
	double ep0 = 1.0 / (exp(p0*h) + 1);
	rhs[0] += (1 - ep0) * f0;
	double pL = eval_func(a, L-h) / eval_func(k, L-h);
	double epL = 1.0 / (exp(pL*h) + 1);
	rhs[n-1] += epL * fL;

	f_out[0] = f0;
	f_out[N-1] = fL;
	nb_matrix_solve_tridiagonal_destructive(n, ld, md, ud, rhs,
						&(f_out[1]));

	nb_soft_free_mem(memsize, memblock);

EXIT:
	return;
}

void nb_analytic_fem_solve_diffusion_convection
			(uint32_t N, double L,
			 double (*k)(double x), /* Diffusion (can be NULL) */
			 double (*a)(double x), /* Convection (can be NULL) */
			 double (*c)(double x), /* Source (can be NULL) */
			 double f0, double fL,  /* Boundary conditions */
			 double *f_out)
/* -k(x) f'' + a(x) f' = c(x) */
{
	if (N < 3)
		goto EXIT;

	uint32_t n = N - 2;
	uint32_t memsize = 4 * n * sizeof(double);
	char *memblock = nb_soft_allocate_mem(memsize);
	double *ld = (void*) memblock;
	double *md = (void*) (memblock + n * sizeof(double));
	double *ud = (void*) (memblock + 2 * n * sizeof(double));
	double *rhs = (void*) (memblock + 3 * n * sizeof(double));

	double h = L/(N-1);
	for (uint32_t i = 0; i < n; i++) {
		double x = h + i * h;
		double kx = eval_func(k, x);
		double ax = eval_func(a, x);
		double cx = eval_func(c, x);
		double p = ax / kx;
		double aux = exp(p*h) - 1;
		double ep = 1.0 / (exp(p*h) + 1);
		ld[i] = - (1 - ep);
		md[i] = 1.0;
		ud[i] = - (ep);
		rhs[i] = h * (cx/ax) * (1 - 2*ep);
	}
	double p0 = eval_func(a, h) / eval_func(k, h);
	double ep0 = 1.0 / (exp(p0*h) + 1);
	rhs[0] += (1 - ep0) * f0;
	double pL = eval_func(a, L-h) / eval_func(k, L-h);
	double epL = 1.0 / (exp(pL*h) + 1);
	rhs[n-1] += epL * fL;

	f_out[0] = f0;
	f_out[N-1] = fL;
	nb_matrix_solve_tridiagonal_destructive(n, ld, md, ud, rhs,
						&(f_out[1]));

	nb_soft_free_mem(memsize, memblock);

EXIT:
	return;
}
