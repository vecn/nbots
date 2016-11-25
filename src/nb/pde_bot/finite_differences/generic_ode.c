#include <stdlib.h>
#include <stdint.h>

#include "nb/memory_bot.h"

#include "nb/pde_bot/finite_differences/generic_ode.h"

#define POW2(a) ((a)*(a))

double eval_func(double (*f)(double x), double x);

void nb_fd_solve_ode_h2(uint32_t N, double L,
			double (*k)(double x), /* Diffusion (can be NULL) */
			double (*a)(double x), /* Convection (can be NULL) */
			double (*b)(double x), /* Absorption (can be NULL) */
			double (*c)(double x), /* Source (can be NULL) */
			double f0, double fL,  /* Boundary conditions */
			double *f_out)
/* -k f'' + a f' + b f = c */
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
