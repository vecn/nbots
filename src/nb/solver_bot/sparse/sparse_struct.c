#include <stdlib.h>
#include <stdint.h>

#include "nb/memory_bot.h"

#include "sparse_struct.h"

uint32_t nb_sparse_bsearch_row(const nb_sparse_t *const A,
			       uint32_t i, uint32_t col, 
			       int imin, int imax)
{
	uint32_t index = (imin + imax) / 2;
	while (A->rows_index[i][index] != col && imin <= imax) {
		index = (imin + imax) / 2;
		if (A->rows_index[i][index] < col)
			imin = index + 1;
		else
			imax = index - 1;
	}
	return index;
}

nb_sparse_t* sparse_allocate(uint32_t N)
{
	nb_sparse_t* A = nb_allocate_mem(sizeof(*A));
	A->rows_values = nb_allocate_mem(N * sizeof(*(A->rows_values)));
	A->rows_index = nb_allocate_mem(N * sizeof(*(A->rows_index)));
	A->rows_size = nb_allocate_zero_mem(N * sizeof(*(A->rows_size)));
	A->N = N;
	return A;
}
