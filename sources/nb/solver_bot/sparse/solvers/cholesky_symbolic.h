#ifndef __NB_SOLVER_BOT_SPARSE_SOLVERS_CHOLESKY_SYMBOLIC_H__
#define __NB_SOLVER_BOT_SPARSE_SOLVERS_CHOLESKY_SYMBOLIC_H__

#include "nb/solver_bot/sparse/sparse.h"

void nb_sparse_cholesky_symbolic(const nb_sparse_t *const A,
				 nb_sparse_t *_L, nb_sparse_t* _Lt, 
				 uint32_t ktrunc);

#endif
