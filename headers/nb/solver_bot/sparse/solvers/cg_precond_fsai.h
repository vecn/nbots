#ifndef __NB_SOLVER_BOT_SPARSE_SOLVERS_CG_PRECOND_FSAI_H__
#define __NB_SOLVER_BOT_SPARSE_SOLVERS_CG_PRECOND_FSAI_H__

#include <stdint.h>

#include "nb/solver_bot/sparse/sparse.h"

int nb_sparse_solve_CG_precond_fsai
(const nb_sparse_t *const A,
 const double *const b,
 double *_x,                /* Out */
 double threshold,
 uint32_t max_iter,	double tolerance,
 uint32_t* niter_performed,     /* Out (NULL if not required) */
 double* tolerance_reached, /* Out (NULL if not required) */
 uint32_t omp_parallel_threads);

#endif
