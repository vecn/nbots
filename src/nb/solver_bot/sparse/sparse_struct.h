#ifndef __NB_EIGEN_BOT_SPARSE_SPARSE_STRUCT_H__
#define __NB_EIGEN_BOT_SPARSE_SPARSE_STRUCT_H__

struct vcn_sparse_s {
	double **rows_values;
	uint32_t **rows_index;
	uint32_t *rows_size;
	uint32_t N;
};

#endif
