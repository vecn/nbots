#ifndef __NB_SOLVER_BOT_SPARSE_EIGEN_LANCZOS_H__
#define __NB_SOLVER_BOT_SPARSE_EIGEN_LANCZOS_H__

#include <stdint.h>

#include "nb/solver_bot/sparse/sparse.h"

void nb_sparse_eigen_lanczos(const nb_sparse_t *const A,
			      double *_eigenmax,/* Out */ 
			      double* _eigenmin,/* Out */
			      int* it,          /* Out */
			      double tolerance,
			      uint32_t omp_parallel_threads);

#endif
