#include <stdint.h>
#include <math.h>

#include "nb/solver_bot/vector.h"

#define POW2(a) ((a)*(a))

double nb_vector_get_norm(const double* x, uint32_t N)
{
	double n = 0;
	for (uint32_t i = 0; i < N; i++)
		n += POW2(x[i]);
	return sqrt(n);
}

void nb_vector_permutation(uint32_t N, const double *v,
			   const uint32_t *perm, double *vp)
{
	for (uint32_t i = 0; i < N; i++)
		vp[i] = v[perm[i]];
}

void nb_vector_sum(uint32_t N, double *a, const double *b)
{
	for (uint32_t i = 0; i < N; i++)
		a[i] += b[i];
}

void nb_vector_substract(uint32_t N, double *a, const double *b)
{
	for (uint32_t i = 0; i < N; i++)
		a[i] -= b[i];
}

void nb_vector_substract_to(uint32_t N, double *a, const double *b)
{
	for (uint32_t i = 0; i < N; i++)
		a[i] = b[i] - a[i];
}

void nb_vector_scale(uint32_t N, double *a, double factor)
{
	for (uint32_t i = 0; i < N; i++)
		a[i] *= factor;
}
