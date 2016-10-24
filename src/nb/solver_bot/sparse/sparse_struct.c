#include <stdlib.h>
#include <stdint.h>

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
