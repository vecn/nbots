#ifndef __NB_SOLVER_BOT_SPARSE_EIGEN_POWER_H__
#define __NB_SOLVER_BOT_SPARSE_EIGEN_POWER_H__

#include <stdint.h>

#include "nb/solver_bot/sparse/sparse.h"


void nb_sparse_eigen_power(const nb_sparse_t *const A, int h,
			   double **_eigenvecs,/* Out */
			   double* _eigenvals, /* Out */
			   int* it,            /* Out */
			   double tolerance,
			   uint32_t omp_parallel_threads);

#endif
