#ifndef __NB_SOLVER_BOT_H__
#define __NB_SOLVER_BOT_H__

#include "nb/solver_bot/vector.h"
#include "nb/solver_bot/matrix/matrix2X2.h"
#include "nb/solver_bot/matrix/matrix3X3.h"
#include "nb/solver_bot/matrix/triangular.h"
#include "nb/solver_bot/matrix/cholesky.h"
#include "nb/solver_bot/matrix/cond.h"
#include "nb/solver_bot/matrix/svd.h"
#include "nb/solver_bot/matrix/qr.h"
#include "nb/solver_bot/sparse/sparse.h"
#include "nb/solver_bot/sparse/spyplot.h"
#include "nb/solver_bot/sparse/solvers/gauss_seidel.h"
#include "nb/solver_bot/sparse/solvers/conjugate_gradient.h"
#include "nb/solver_bot/sparse/solvers/cg_precond_jacobi.h"
#include "nb/solver_bot/sparse/solvers/cg_precond_chol.h"
#include "nb/solver_bot/sparse/solvers/cg_precond_fsai.h"
#include "nb/solver_bot/sparse/solvers/triangular.h"
#include "nb/solver_bot/sparse/solvers/cholesky.h"
#include "nb/solver_bot/sparse/solvers/lu.h"
#include "nb/solver_bot/sparse/eigen/power.h"
#include "nb/solver_bot/sparse/eigen/inv_power.h"
#include "nb/solver_bot/sparse/eigen/lanczos.h"
#include "nb/solver_bot/sparse/eigen/givens.h"
#include "nb/solver_bot/matlab_v4.h"

#endif
