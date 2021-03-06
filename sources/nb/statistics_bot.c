#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

#include "nb/memory_bot.h"
#include "nb/math_bot.h"
#include "nb/container_bot/array.h"
#include "nb/statistics_bot.h"

inline uint32_t nb_statistics_get_seed(void)
{
	srand(time(0));
	return rand();	
}

inline uint32_t nb_statistics_lcg(uint32_t seed)
/* Linear congruencial generator (Park and Miller 1988) */
{
	if (0 == seed)
		seed = 4294967291UL; /* 2^32 - 5 */
	return ((uint64_t)seed * 48271UL) % 2147483647UL;
}

void nb_statistics_random_permutation(uint32_t N, void *base, 
				      uint16_t type_size)
{	
	uint32_t rseed = nb_statistics_get_seed();
	for (uint32_t i = 0; i < N; i++) {
		uint32_t k = rseed % N;
		rseed = nb_statistics_lcg(rseed);
		nb_swap(base, i, k, type_size);
	}
}

void nb_statistics_runif(int n, double min, double max, 
			  double *const restrict out,
			  uint64_t *const restrict seed)
{
	/* Generate n random numbers from a uniform distribution
	 * using the Park and Miller (1988) algorithm.
	 */
	/* A = 16807 */
	/* M = 2^31 -1 = 2147483647 */
	/* Q = M / A = 127773 */
	/* R = M % A = 2836 */
	for (uint32_t i = 0; i < n; i++) {
		/*(A)*/           /*(Q)*/  /*(R)*/          /*(Q)*/
		*seed = 16807 * (*seed % 127773) - 2836 * (*seed / 127773);
		/* (M-1) */
		out[i] = min + (max-min) * (*seed/2147483647.0);
	}
}

void nb_statistics_rnorm(int n, double mean, double var,
			  double *const restrict out,
			  uint64_t *const restrict seed1, 
			  uint64_t *const restrict seed2)
{
	/* Generate n random numbers from a normal distribution
	 * using the Box-Muller transformation method.
	 */
	const int m = (int)(n * 0.5 + 0.5);
	double *const restrict U = nb_allocate_mem(m * sizeof(double));
	double *const restrict V = nb_allocate_mem(m * sizeof(double));
	nb_statistics_runif(m, 0.0, 1.0, U, seed1);
	nb_statistics_runif(m, 0.0, 1.0, V, seed2);
	for (uint32_t i = 0; i < m; i++) {
		out[i*2] = mean + sqrt(-2.0 * log(U[i])) * 
			cos(2.0 * NB_MATH_PI * V[i]) * var;
		if (i*2+1 < n) 
			out[i*2+1] = mean + sqrt(-2 * log(V[i])) * 
				sin(2.0 * NB_MATH_PI * U[i]) * var;
	}
	nb_free_mem(U);
	nb_free_mem(V);
}
