#ifndef __NB_SOLVER_BOT_SPARSE_SOLVERS_TRIANGULAR_H__
#define __NB_SOLVER_BOT_SPARSE_SOLVERS_TRIANGULAR_H__

#include <stdint.h>

#include "nb/solver_bot/sparse/sparse.h"


void nb_sparse_solve_LU(const nb_sparse_t *const L, 
			 const nb_sparse_t *const U,
			 const double *const b,
			 double* _x);  /* Out */

void nb_sparse_forward_solve(const nb_sparse_t *const L,
			      const double *const b, 
			      double* _x /* Out */);
void nb_sparse_backward_solve(const nb_sparse_t *const U,
			       const double *const b,
			       double* _x /* Out */);

double nb_sparse_triangular_get_det_sign(const nb_sparse_t *Trg);

#endif
