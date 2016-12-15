#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include "nb/memory_bot.h"
#include "nb/solver_bot.h"

#define POW2(a) ((a)*(a))

static void calculate_G_and_sum_Ah(uint32_t N, uint8_t dim_x, uint8_t dim_f,
				   const double *x, const double *f,
				   const double *ni, const double *fi,
				   double *G, double *Ah_sum);
static double get_normal(uint8_t dim_x, const double *ni,
			 uint32_t id, const double *x, double *nij);
static void get_fij(uint8_t dim_f, const double *fi, const double *f,
		    uint32_t i, double *fij);

static void sum_G(uint8_t dim_x, const double *nij, double *G);
static void sum_Ah(uint8_t dim_x, uint8_t dim_f,
		   const double *nij, const double *fij,
		   double dist, double *Ah_sum);
static void transpose_Df(double *Df, double *Df_aux,
			 uint8_t dim_x, uint8_t dim_f);

void nb_pde_get_frechet_derivative(uint32_t N, uint8_t dim_x, uint8_t dim_f,
				   const double *x, const double *f,
				   const double *ni, const double *fi,
				   double *Df)
{
	uint16_t G_size = POW2(dim_x) * sizeof(double);
	uint16_t Ah_size = dim_x * dim_f * sizeof(double);
	uint32_t memsize = 2 * G_size + Ah_size;
	char *memblock = nb_soft_allocate_mem(memsize);
	double *G = (void*) memblock;
	double *LLt = (void*) (memblock + G_size);
	double *Ah_sum = (void*) (memblock + 2 * G_size);

	calculate_G_and_sum_Ah(N, dim_x, dim_f, x, f, ni, fi, G, Ah_sum);

	int status = nb_matrix_cholesky_decomposition(G, LLt, dim_x);
  
	if (0 != status) {
		/* Single component information */
		for (uint16_t i = 0; i < dim_x * dim_f; i++)
			Df[i] = Ah_sum[i] / N;
	} else {
		for (uint16_t i = 0; i < dim_f; i++)
			nb_matrix_cholesky_solve(LLt, &(Ah_sum[i*dim_x]),
						  &(Df[i*dim_x]), dim_x);
	}
	transpose_Df(Df, Ah_sum, dim_x, dim_f);
	nb_soft_free_mem(memsize, memblock);
}

static void calculate_G_and_sum_Ah(uint32_t N, uint8_t dim_x, uint8_t dim_f,
				   const double *x, const double *f,
				   const double *ni, const double *fi,
				   double *G, double *Ah_sum)
{
	uint16_t memsize = (dim_x + dim_f) * sizeof(double);
	char *memblock = nb_soft_allocate_mem(memsize);
	double *nij = (void*) memblock;
	double *fij = (void*) (memblock + dim_x * sizeof(double));

	memset(G, 0, POW2(dim_x) * sizeof(*G));
	memset(Ah_sum, 0, dim_x * dim_f * sizeof(*Ah_sum));
	for (uint32_t i = 0; i < N; i++) {
		double dist = get_normal(dim_x, ni, i, x, nij);
		get_fij(dim_f, fi, f, i, fij);
		sum_G(dim_x, nij, G);
		sum_Ah(dim_x, dim_f, nij, fij, dist, Ah_sum);
	}
	nb_soft_free_mem(memsize, memblock);
}

static double get_normal(uint8_t dim_x, const double *ni,
			 uint32_t id, const double *x, double *nij)
{
	double dist = 0;
	for (uint8_t d = 0; d < dim_x; d++) {
		nij[d] = ni[id * dim_x + d] - x[d];
		dist += POW2(nij[d]);
	}
	dist = sqrt(dist);
	for (uint8_t d = 0; d < dim_x; d++)
		nij[d] /= dist;
	return dist;
}

static void get_fij(uint8_t dim_f, const double *fi, const double *f,
		    uint32_t i, double *fij)
{
	for (uint8_t d = 0; d < dim_f; d++)
		fij[d] = fi[i * dim_f + d] - f[d];
}

static void sum_G(uint8_t dim_x, const double *nij, double *G)
{
	for (uint8_t l = 0; l < dim_x; l++) {
		for (uint8_t m = 0; m < dim_x; m++)
			G[l * dim_x + m] += nij[l] * nij[m];
	}
}

static void sum_Ah(uint8_t dim_x, uint8_t dim_f,
		   const double *nij, const double *fij,
		   double dist, double *Ah_sum)
{
	for (uint8_t l = 0; l < dim_f; l++) {
		for (uint8_t m = 0; m < dim_x; m++)
			Ah_sum[l * dim_x + m] += (fij[l] * nij[m]) / dist;
	}
}

static void transpose_Df(double *Df, double *Df_aux,
			 uint8_t dim_x, uint8_t dim_f)
{
	memcpy(Df_aux, Df, dim_x * dim_f * sizeof(double));
	for (uint8_t i = 0; i < dim_x; i++) {
		for (uint8_t j = 0; j < dim_f; j++) {
			Df[i * dim_f + j] = Df_aux[j * dim_x + i];
		}
	}
}
