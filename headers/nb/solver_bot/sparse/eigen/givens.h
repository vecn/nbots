#ifndef __NB_SOLVER_BOT_SPARSE_EIGEN_GIVENS_H__
#define __NB_SOLVER_BOT_SPARSE_EIGEN_GIVENS_H__

#include <stdbool.h>
#include <stdint.h>

#include "nb/solver_bot/sparse/sparse.h"

void nb_sparse_eigen_givens(const double* const main_diag, 
			     const double* const uplw_diag,
			     int i, double *_eigenvalue,
			     double tolerance,
			     uint32_t N);
#endif
