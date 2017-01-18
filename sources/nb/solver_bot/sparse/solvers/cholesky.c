#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "nb/math_bot.h"
#include "nb/memory_bot.h"
#include "nb/container_bot.h"
#include "nb/graph_bot.h"
#include "nb/solver_bot.h"

#include "../sparse_struct.h"

#define POW2(a) ((a)*(a))

int nb_sparse_decompose_Cholesky(const nb_sparse_t *const A,
				  nb_sparse_t *L,             /* Out */
				  nb_sparse_t* Lt,            /* Out */
				  uint32_t omp_parallel_threads)
{
	/* "L" must be a lower triangular matrix with the main diagonal 
	 * complete, and "Lt" must be an upper triangular matrix with
	 * the main diagonal complete.
	 * The structure of L must be congrous with Lt, since Lt = L'.
	 */
    
	/* Compute the decomposition */
	for(uint32_t j=0; j< A->N; j++){
		L->rows_values[j][L->rows_size[j]-1] = nb_sparse_get(A, j, j);
    
		double sum = 0;
		for(uint32_t q = 0; q < L->rows_size[j]-1; q++)
			sum += POW2(L->rows_values[j][q]);
      
		L->rows_values[j][L->rows_size[j]-1] -= sum;

		if(L->rows_values[j][L->rows_size[j]-1] <= 0.0)
			return 1;
              
		double valuejj = sqrt(L->rows_values[j][L->rows_size[j]-1]);
		L->rows_values[j][L->rows_size[j]-1] = valuejj;
		Lt->rows_values[j][0] = valuejj;

#pragma omp parallel for  num_threads(omp_parallel_threads) schedule(guided)
		for(uint32_t q = 1; q < Lt->rows_size[j]; q++){
			uint32_t i = Lt->rows_index[j][q];
			/*** L_ij <- A_ij ********************************************************/
			uint32_t L_jindex = nb_sparse_bsearch_row(L, i, j, 0, L->rows_size[i]-1);     /**/
			L->rows_values[i][L_jindex] = nb_sparse_get(A, i, j);                     /**/
			/*************************************************************************/
			uint32_t r = 0;
			uint32_t s = 0;
			uint32_t _ro = L->rows_index[i][r];
			uint32_t _sigma = L->rows_index[j][s];
			bool flag = true;   /* Flag to know when to stop the cylce */
			while(flag){
				while(_ro < _sigma)
					_ro = L->rows_index[i][++r];
				while(_ro > _sigma)
					_sigma = L->rows_index[j][++s];
				while(_ro == _sigma){
					if(_ro == j){
						flag = false;   /* Finish the cycle */
						break;
					}
					double vir = L->rows_values[i][r];
					double vjs = L->rows_values[j][s];
					L->rows_values[i][L_jindex] -= vir*vjs;
					_ro = L->rows_index[i][++r];
					_sigma = L->rows_index[j][++s];
				}
			}

			L->rows_values[i][L_jindex] /= L->rows_values[j][L->rows_size[j]-1];
			Lt->rows_values[j][q] = L->rows_values[i][L_jindex];
		}
	}

	/* Successful exit */
	return 0;
}

int nb_sparse_solve_Cholesky(const nb_sparse_t *const A,
			      const double *const b,
			      double* x,  /* Out */
			      uint32_t omp_parallel_threads){
	nb_sparse_t *L = NULL; 
	nb_sparse_t *U = NULL;
	nb_sparse_alloc_LU(A, &L, &U);
	if (NULL == L)
		return 10;
	
	int status = nb_sparse_decompose_Cholesky(A, L, U, 
						  omp_parallel_threads);

	if (0 == status)
		nb_sparse_solve_LU(L, U, b, x);

	nb_sparse_destroy(L);
	nb_sparse_destroy(U);

	return status;
}

int nb_sparse_relabel_and_solve_Cholesky(const nb_sparse_t *const A,
					       const double *const b,
					       double* x,  /* Out */
					       uint32_t omp_parallel_threads)
{
	uint32_t N = nb_sparse_get_size(A);
	uint32_t memsize = 2 * N * (sizeof(uint32_t) + sizeof(double));
	char *memblock = nb_soft_allocate_mem(memsize);
	uint32_t *perm = (void*) memblock;
	uint32_t *iperm = (void*) (memblock + N * sizeof(uint32_t));
	double *br = (void*) (memblock + 2 * N * sizeof(uint32_t));
	double *xr = (void*) (memblock + 2 * N * sizeof(uint32_t) +
			      N * sizeof(double));

	nb_sparse_calculate_permutation(A, perm, iperm);

	nb_sparse_t *Ar = nb_sparse_create_permutation(A, perm, iperm);
	nb_vector_permutation(N, b, perm, br);

	int status = nb_sparse_solve_Cholesky(Ar, br, xr,
					      omp_parallel_threads);

	nb_vector_permutation(N, xr, iperm, x);
	
	nb_sparse_destroy(Ar);
	nb_soft_free_mem(memsize, memblock);
	return status;
}
