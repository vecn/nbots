#ifndef __NB_SOLVER_BOT_SPARSE_EIGEN_INV_POWER_H__
#define __NB_SOLVER_BOT_SPARSE_EIGEN_INV_POWER_H__

#include <stdint.h>

#include "nb/solver_bot/sparse/sparse.h"

int nb_sparse_eigen_ipower(const nb_sparse_t *const A,
			    nb_solver_t solver,
			    int h, double mu,
			    double **_eigenvecs,/* Out */
			    double* _eigenvals, /* Out */
			    int* it,            /* Out */
			    double tolerance,
			    uint32_t omp_parallel_threads);

#endif
