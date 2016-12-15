#include <stdint.h>

#include "nb/memory_bot.h"

#include "nb/solver_bot/matrix/tridiagonal_solver.h"

void nb_matrix_solve_tridiagonal_destructive(uint32_t N, const double *ld,
					     const double *md, double *ud,
					     const double *b, double *x)
{
	ud[0] = ud[0] / md[0];
	x[0] = b[0] / md[0];

	for (uint32_t i = 1; i < N; i++) {
		double m = 1.0 / (md[i] - ld[i] * ud[i - 1]);
		ud[i] = ud[i] * m;
		x[i] = (b[i] - ld[i] * x[i - 1]) * m;
	}
	for (int32_t i = N - 2; i >= 0; i--)
		x[i] = x[i] - ud[i] * x[i + 1];
}

void nb_matrix_solve_tridiagonal(uint32_t N, const double *ld,
				 const double *md, const double *ud,
				 const double *b, double *x)
{
	uint32_t memsize = N * sizeof(double);
	double *aux = nb_soft_allocate_mem(memsize);

	aux[0] = ud[0] / md[0];
	x[0] = b[0] / md[0];

    	for (uint32_t i = 1; i < N; i++) {
		double m = 1.0 / (md[i] - ld[i] * aux[i - 1]);
		aux[i] = ud[i] * m;
		x[i] = (b[i] - ld[i] * x[i - 1]) * m;
	}

    	for (int32_t i = N - 2; i >= 0; i--)
		x[i] = x[i] - aux[i] * x[i + 1];

    	nb_soft_free_mem(memsize, aux);
}
