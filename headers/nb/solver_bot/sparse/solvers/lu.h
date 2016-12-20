#ifndef __NB_SOLVER_BOT_SPARSE_SOLVERS_LU_H__
#define __NB_SOLVER_BOT_SPARSE_SOLVERS_LU_H__

#include <stdbool.h>
#include <stdint.h>

#include "nb/solver_bot/sparse/sparse.h"

/**
 * @brief LU decomposition (Doolitle).
 * L and U must be allocated previously using nb_sparse_alloc_LU().
 * @param[out] L Lower triangular matrix resulting from the decomposition.
 * @param[out] U Upper triangular matrix resulting from the decomposition.
 * @param[in] Number of threads to be used in parallel task.
 */
void nb_sparse_decompose_LU(const nb_sparse_t *const Ar,
			     nb_sparse_t *L, nb_sparse_t* U,
			     uint32_t omp_parallel_threads);


int nb_sparse_solve_using_LU(const nb_sparse_t *const A,
			      const double *const b,
			      double* x,  /* Out */
			      uint32_t omp_parallel_threads);

int nb_sparse_relabel_and_solve_using_LU(const nb_sparse_t *const A,
					 const double *const b,
					 double* x,  /* Out */
					 uint32_t omp_parallel_threads);

#endif
