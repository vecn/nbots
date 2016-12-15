#ifndef __NB_SOLVER_BOT_SPARSE_SOLVERS_CHOLESKY_H__
#define __NB_SOLVER_BOT_SPARSE_SOLVERS_CHOLESKY_H__

#include <stdint.h>

#include "nb/solver_bot/sparse/sparse.h"

/**
 * @brief Allocate LU decomposition using symbolic Cholesky
 * @param[in] A Matrix to be decomposed.
 * @param[in] L Pointer to the Lower triangular matrix that will be allocated.
 * @param[in] U Pointer to the Upper triangular matrix that will be allocated.
 * @return Zero if success. The matrix will be initialized in zeros.
 */
int nb_sparse_alloc_LU
(const nb_sparse_t *const A, nb_sparse_t** L, nb_sparse_t** U);

/**
 * @brief Cholesky decomposition.
 * L and Lt must be allocated previously using nb_sparse_alloc_LU().
 * @param[out] L Lower triangular matrix resulting from the decomposition.
 * @param[out] Lt Upper triangular matrix resulting from the decomposition.
 * @param[in] Number of threads to be used in parallel task.
 */
int nb_sparse_decompose_Cholesky(const nb_sparse_t *const Ar, 
				  nb_sparse_t * L, nb_sparse_t* Lt,
				  uint32_t omp_parallel_threads);
  
int nb_sparse_solve_Cholesky(const nb_sparse_t *const A,
			      const double *const b,
			      double* x,  /* Out */
			      uint32_t omp_parallel_threads);


#endif
