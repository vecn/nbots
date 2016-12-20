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
#include "cholesky_symbolic.h"

int nb_sparse_alloc_LU(const nb_sparse_t *const restrict A,
		       nb_sparse_t** L, nb_sparse_t** U)
{
	*L = nb_sparse_allocate(A->N);
	*U = nb_sparse_allocate(A->N);
  
	nb_sparse_cholesky_symbolic(A, *L, *U, A->N);

	return 0;
}

void nb_sparse_decompose_LU(const nb_sparse_t *const Ar,
			     nb_sparse_t *L, nb_sparse_t* U,
			     uint32_t omp_parallel_threads)
{

	/* Create Ut to compute faster the decomposition */
	nb_sparse_t* Ut = nb_sparse_clone(L);

	/* Compute the decomposition */
	for(uint32_t j=0; j< Ar->N; j++){
		L->rows_values[j][L->rows_size[j]-1] = 1.0;
		U->rows_values[j][0] = nb_sparse_get(Ar, j, j);

		double sum = 0;

#pragma omp parallel for schedule(guided) reduction(+:sum) num_threads(omp_parallel_threads)
		for(uint32_t q=0; q< L->rows_size[j]-1; q++)
			sum += L->rows_values[j][q] * Ut->rows_values[j][q];
    

		U->rows_values[j][0] -= sum;
		Ut->rows_values[j][Ut->rows_size[j]-1] = U->rows_values[j][0];
        
#pragma omp parallel for schedule(guided) num_threads(omp_parallel_threads)
		for(uint32_t q = 1; q < U->rows_size[j]; q++){
			uint32_t i = U->rows_index[j][q];
			/*** L_ij <- A_ij *******************************************************/
			uint32_t L_jindex = nb_sparse_bsearch_row(L, i, j, 0, L->rows_size[i]-1);/**/
			L->rows_values[i][L_jindex] = nb_sparse_get(Ar, i, j);               /**/
			/************************************************************************/
			/*** U_ji <- A_ji *******************************************************/
			U->rows_values[j][q] = nb_sparse_get(Ar, j, i);                      /**/
			/************************************************************************/
			register uint32_t r = 0;
			register uint32_t s = 0;
			register uint32_t _ro =    L->rows_index[i][r];
			register uint32_t _sigma = L->rows_index[j][s];
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
					double vjs = Ut->rows_values[j][s];
					L->rows_values[i][L_jindex] -= vir*vjs;

					vjs = L->rows_values[j][s];
					vir = Ut->rows_values[i][r];
					U->rows_values[j][q] -= vir*vjs;

					_ro = L->rows_index[i][++r];
					_sigma = L->rows_index[j][++s];
				}
			}

			L->rows_values[i][L_jindex] /= U->rows_values[j][0];
      
			Ut->rows_values[i][L_jindex] = U->rows_values[j][q];
		}
	}
	/* Free memory */
	nb_sparse_destroy(Ut);
}

void  nb_sparse_solve_LU(const nb_sparse_t *const L, 
			  const nb_sparse_t *const U,
			  const double *const b,
			  double* _x  /* Out */)
{
	double* z = nb_allocate_zero_mem(L->N * sizeof(*z));
	nb_sparse_forward_solve(L, b, z);
	nb_sparse_backward_solve(U, z, _x);
	nb_free_mem(z);
}

int nb_sparse_solve_using_LU(const nb_sparse_t *const A,
			      const double *const b,
			      double* x,  /* Out */
			      uint32_t omp_parallel_threads)
{
	nb_sparse_t *L = NULL; 
	nb_sparse_t *U = NULL;
	nb_sparse_alloc_LU(A, &L, &U);
	if(NULL == L)
		return 1;

	nb_sparse_decompose_LU(A, L, U, omp_parallel_threads);

	nb_sparse_solve_LU(L, U, b, x);

	nb_sparse_destroy(L);
	nb_sparse_destroy(U);

	return 0;
}

int nb_sparse_relabel_and_solve_using_LU(const nb_sparse_t *const A,
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

	int status = nb_sparse_solve_using_LU(Ar, br, xr,
					      omp_parallel_threads);

	nb_vector_permutation(N, xr, iperm, x);
	
	nb_sparse_destroy(Ar);
	nb_soft_free_mem(memsize, memblock);
	return status;
}
