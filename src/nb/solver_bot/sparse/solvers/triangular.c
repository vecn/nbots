#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "nb/solver_bot/sparse/sparse.h"
#include "nb/solver_bot/sparse/solvers/triangular.h"

#include "../sparse_struct.h"

void nb_sparse_forward_solve(const nb_sparse_t *const L,
			      const double *const b, 
			      double* _x /* Out */)
/* Solve the system Lx = b, where L is a lower triangular matrix.
 * If L is not a lower triangular matrix, the program are going to
 * take only the diagonal and the diagonal-lower part of the matrix.
 */
{
	/* This solver cannot be parallelized to guaranty the solution */
	for(int i=0; i< L->N; i++){
		_x[i] = b[i];
		register int q = 0;
		while(L->rows_index[i][q] < i){
			_x[i] -= L->rows_values[i][q] * _x[L->rows_index[i][q]];
			q++;
		}    
		if(L->rows_index[i][L->rows_size[i]-1] == i)
			_x[i] /= L->rows_values[i][L->rows_size[i]-1];
		else
			_x[i] /= nb_sparse_get(L, i, i);
	}
}

void nb_sparse_backward_solve(const nb_sparse_t *const U,
			       const double *const b,
			       double* _x /* Out */)
/* Solve the system Ux = b, where U is a upper triangular matrix.
 * If U is not an upper triangular matrix, the program are going to
 * take only the diagonal and the diagonal-upper part of the matrix.
 */
{
	/* This solver cannot be parallelized to guaranty the solution */
	for(int i = U->N-1; i >= 0 ; i--){
		_x[i] = b[i];
		uint32_t idx = 0;
		if(U->rows_index[i][0] != i)
			idx = nb_sparse_bsearch_row(U, i, i, 0, U->rows_size[i]-1);
		for(int q = idx+1; q < U->rows_size[i]; q++)
			_x[i] -= U->rows_values[i][q]*_x[U->rows_index[i][q]];
    
		_x[i] /= U->rows_values[i][idx];
	}
}
