#ifndef __NB_SOLVER_BOT_SPARSE_SPARSE_STRUCT_H__
#define __NB_SOLVER_BOT_SPARSE_SPARSE_STRUCT_H__

#include "nb/solver_bot/sparse/sparse.h"

struct nb_sparse_s {
	double **rows_values;
	uint32_t **rows_index;
	uint32_t *rows_size;
	uint32_t N;
};

nb_sparse_t* nb_sparse_allocate(uint32_t N);

uint32_t nb_sparse_bsearch_row(const nb_sparse_t *const A,
			       uint32_t i, uint32_t col, 
			       int imin, int imax);

#endif
