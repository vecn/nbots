#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include "nb/math_bot.h"
#include "nb/memory_bot.h"
#include "nb/container_bot.h"
#include "nb/solver_bot/sparse/sparse.h"

#include "sparse_struct.h"

#define POW2(a) ((a)*(a))

static void verify_graph(const nb_graph_t *graph);
static int meta_compare_data_bycol(const void *a, const void *b);

nb_sparse_t* nb_sparse_create(const nb_graph_t *const restrict graph,
				const uint32_t *const restrict perm,
				uint32_t vars_per_node)
{
	verify_graph(graph);
	nb_sparse_t*A = sparse_allocate(graph->N * vars_per_node);

	for (uint32_t i = 0; i < graph->N; i++) {
		uint32_t row_size = (graph->N_adj[i] + 1) * vars_per_node;    
		for (uint32_t k1 = 0; k1 < vars_per_node; k1++) {
			uint32_t irow;
			if (NULL != perm)
				irow = perm[i] * vars_per_node + k1;
			else 
				irow = i * vars_per_node + k1;

			A->rows_size[irow] = row_size;
			A->rows_values[irow] = 
				nb_allocate_zero_mem(row_size *
				       sizeof(*(A->rows_values[irow])));
			A->rows_index[irow] = 
				nb_allocate_mem(row_size *
				       sizeof(*(A->rows_index[irow])));

			for (uint32_t j = 0; j < graph->N_adj[i]; j++) {
				uint32_t id = graph->adj[i][j];
				for (uint32_t k2 = 0; k2 < vars_per_node; k2++) {
					uint32_t icol = id * vars_per_node + k2;
					A->rows_index[irow][j * vars_per_node + k2] = icol;
				}
			}
			for (uint32_t k2 = 0; k2 < vars_per_node; k2++) {
				uint32_t icol = i * vars_per_node + k2;
				A->rows_index[irow][graph->N_adj[i] * vars_per_node + k2] = icol ;
			}
			nb_qsort(A->rows_index[irow], A->rows_size[irow], 
				  sizeof(uint32_t), nb_compare_uint32);
		}
	}
	return A;
}

static void verify_graph(const nb_graph_t *graph)
{
#ifndef NDEBUG
	uint32_t N = graph->N;
	for (uint32_t i = 0; i < N; i++) {
		for (uint32_t j = 0; j < graph->N_adj[i]; j++) {
			uint32_t adj = graph->adj[i][j];
			assert(adj < N);
			assert(i != adj);
			for (uint32_t k = j + 1; k < graph->N_adj[i]; k++)
				assert(adj != graph->adj[i][k]);
		}
	}
#endif
}

nb_sparse_t* nb_sparse_clone(nb_sparse_t* A)
{
	nb_sparse_t* Acopy = sparse_allocate(A->N);
	for (uint32_t i = 0; i < A->N; i++) {
		Acopy->rows_index[i] = nb_allocate_zero_mem(A->rows_size[i] *
					      sizeof(**(Acopy->rows_index)));
		memcpy(Acopy->rows_index[i], A->rows_index[i], 
		       A->rows_size[i] * sizeof(**(Acopy->rows_index)));

		Acopy->rows_values[i] = nb_allocate_zero_mem(A->rows_size[i] *
					       sizeof(**(Acopy->rows_values)));
		memcpy(Acopy->rows_values[i], A->rows_values[i],
		       A->rows_size[i] * sizeof(**(Acopy->rows_values)));

		Acopy->rows_size[i] = A->rows_size[i];
	}
	return Acopy;
}

void nb_sparse_save(const nb_sparse_t *const A, const char* filename)
{
	FILE *fp = fopen(filename, "w");
	if (fp == NULL) {
		printf("Error: Opening the file %s to save the matrix.\n",
		       filename);
		return;
	}
	fprintf(fp, "%i %i\n", A->N, A->N);
	for (uint32_t i=0; i < A->N; i++) {
		for (uint32_t j=0; j < A->rows_size[i]; j++)
			fprintf(fp, "%i %i %e \n", i, A->rows_index[i][j],
				A->rows_values[i][j]);
	}
	fclose(fp);
}

void nb_sparse_destroy(nb_sparse_t* A)
{
	/* Clear all rows */
	for (uint32_t i = 0; i < A->N; i++) {
		nb_free_mem(A->rows_values[i]);
		nb_free_mem(A->rows_index[i]);
	}
	nb_free_mem(A->rows_values);
	nb_free_mem(A->rows_index);
	nb_free_mem(A->rows_size);
	nb_free_mem(A);
}

void nb_sparse_reset(nb_sparse_t *A)
{
	for(uint32_t i = 0; i < A->N; i++)
		memset(A->rows_values[i], 0,
		       A->rows_size[i] * sizeof(*(A->rows_values)));  
}

void nb_sparse_set(nb_sparse_t *A, uint32_t i, uint32_t j, double value)
{
	uint32_t index = nb_sparse_bsearch_row(A, i, j, 0, A->rows_size[i]-1);
	if (A->rows_index[i][index] == j) {
		A->rows_values[i][index] = value;
	} else {
		/* OPPORTUNITY: Send this to some logfile */
		printf("ERROR: Entry of sparse matrix is not allocated.\n");
		printf("    -> nb_sparse_set(*A=%p, i=%i, j=%i, value=%lf)\n",
		       (void*) A, i, j, value);
		exit(1);
	}
}

void nb_sparse_set_identity_row(nb_sparse_t* A, uint32_t row)
{
	memset(A->rows_values[row], 0,
	       A->rows_size[row] * sizeof(**(A->rows_values)));
	nb_sparse_set(A, row, row, 1.0);
}

void nb_sparse_make_diagonal(nb_sparse_t* A, double diag_val)
{
	for (uint32_t i = 0; i < A->N; i++){
		memset(A->rows_values[i], 0, 
		       A->rows_size[i] * sizeof(**(A->rows_values)));
		nb_sparse_set(A, i, i, diag_val);
	}
}

double nb_sparse_get(const nb_sparse_t *const A, uint32_t i, uint32_t j)
{
	uint32_t index = nb_sparse_bsearch_row(A, i, j, 0, A->rows_size[i]-1);
	if (A->rows_index[i][index] == j)
		return A->rows_values[i][index];
	else
		return 0.0;
}

double nb_sparse_get_and_set(nb_sparse_t *A, uint32_t i, uint32_t j, double value)
{
	uint32_t index = nb_sparse_bsearch_row(A, i, j, 0, A->rows_size[i]-1);
	if (A->rows_index[i][index] == j) {
		double val = A->rows_values[i][index];
		A->rows_values[i][index] = value;
		return val;
	} else {
		/* OPPORTUNITY: Send this to some logfile */
		printf("ERROR: Entry of sparse matrix is not allocated.\n");
		printf("    -> nb_sparse_get_and_set(*A=%p, i=%i, j=%i, value=%lf)\n",
		       (void*) A, i, j, value);
		exit(1);
	}
}

bool nb_sparse_is_non_zero(const nb_sparse_t *const A, uint32_t i, uint32_t j)
{
	uint32_t index = nb_sparse_bsearch_row(A, i, j, 0, A->rows_size[i]-1);
	if(A->rows_index[i][index] == j)
		return true;
	else
		return false;
}

uint32_t nb_sparse_memory_used(const nb_sparse_t *const A)
{
	uint32_t size = sizeof(short) + 3 * sizeof(uint32_t) + 3*sizeof(void*);
	size += A->N * (2*sizeof(void*) + sizeof(uint32_t));
	for (uint32_t i = 0; i < A->N; i++)
		size += A->rows_size[i] * (sizeof(uint32_t) + sizeof(double));
	return size;
}

void nb_sparse_add(nb_sparse_t *A, uint32_t i, uint32_t j, double value)
{
	uint32_t index = nb_sparse_bsearch_row(A, i, j, 0, A->rows_size[i]-1);
	if (A->rows_index[i][index] == j) {
		A->rows_values[i][index] += value;
	} else {
		/* OPPORTUNITY: Send this to some logfile */
		printf("ERROR: Entry of sparse matrix is not allocated.\n");
		printf("    -> nb_sparse_add(*A=%p, i=%i, j=%i, value=%lf)\n",
		       (void*) A, i, j, value);
		exit(1);
	}
}

void nb_sparse_scale(nb_sparse_t *A, double factor)
{
	for (uint32_t i = 0; i < A->N; i++) {
		for (uint32_t j = 0; j < A->rows_size[i]; j++)
			A->rows_values[i][j] *= factor;
	}
}
 
nb_sparse_t* nb_sparse_create_permutation
(const nb_sparse_t *const A,
 const uint32_t *const perm,
 const uint32_t *const iperm)
{
	/* Use the vectors perm and iperm to compute Ar = PAP'. */ 
	nb_sparse_t* _Ar = sparse_allocate(A->N);
	for (uint32_t i = 0; i < A->N; i++) {
		uint32_t j = perm[i];
		_Ar->rows_size[i] = A->rows_size[j];
		_Ar->rows_index[i] = 
			nb_allocate_zero_mem(_Ar->rows_size[i] *
					     sizeof(**(_Ar->rows_index)));
		_Ar->rows_values[i] =
			nb_allocate_zero_mem(_Ar->rows_size[i] *
					     sizeof(**(_Ar->rows_values)));
		double** data2sort =
			nb_allocate_mem(_Ar->rows_size[i] * sizeof(*data2sort));
		for (uint32_t k = 0; k < A->rows_size[j]; k++) {
			uint32_t m = iperm[A->rows_index[j][k]];
			data2sort[k] =
				nb_allocate_zero_mem(2 * sizeof(**data2sort));
			data2sort[k][0] = (double)m; /* OPPORTUNITY */
			data2sort[k][1] = A->rows_values[j][k];
		}
		qsort(data2sort, _Ar->rows_size[i], sizeof(*data2sort), 
		      meta_compare_data_bycol);/* TEMPORAL: use nb_qsort*/
		for (uint32_t k=0; k< A->rows_size[j]; k++) {
			_Ar->rows_index[i][k] = (uint32_t)data2sort[k][0];
			_Ar->rows_values[i][k] = data2sort[k][1];
		}
		/* Free memory */
		for (uint32_t k=0; k< A->rows_size[j]; k++)
			nb_free_mem(data2sort[k]);
		nb_free_mem(data2sort);
	}
	return _Ar;
}

static int meta_compare_data_bycol(const void *a, const void *b)
{
	if((*(double**)a)[0] == (*(double**)b)[0])
		return (*(double**)a)[1] - (*(double**)b)[1];
	else
		return (*(double**)a)[0] - (*(double**)b)[0];
}

void nb_sparse_fill_permutation(const nb_sparse_t *const A, 
				 nb_sparse_t* Ar,
				 const uint32_t *const perm,
				 const uint32_t *const iperm)
{
	/* Use the vectors perm and iperm to compute Ar = PAP'. */
	for (uint32_t i=0; i < A->N; i++) {
		uint32_t j = perm[i];
		for (uint32_t k = 0; k < A->rows_size[j]; k++) {
			uint32_t m = iperm[A->rows_index[j][k]];
			uint32_t idx = nb_sparse_bsearch_row(Ar, i, m, 0, Ar->rows_size[i]-1);
			Ar->rows_values[i][idx] = A->rows_values[j][k];
		}
	}
}

double* nb_sparse_create_vector_permutation
(const double *const b, 
 const uint32_t *const perm,
 uint32_t N){
	double* br = (double*)nb_allocate_mem(N * sizeof(double));
	for(uint32_t i=0; i < N; i++)
		br[i] = b[perm[i]];
	return br;
}

void nb_sparse_get_transpose(const nb_sparse_t *A, nb_sparse_t *_At)
/* _At must be allocated and initialized */
{
	uint32_t memsize = A->N * sizeof(uint32_t);
	uint32_t *jcount = nb_soft_allocate_mem(memsize);
	memset(jcount, 0, memsize);
	for (uint32_t i = 0; i < A->N; i++) {
		for (uint32_t q = 0; q < A->rows_size[i]; q++) {
			uint32_t j = A->rows_index[i][q];
			uint32_t jc = jcount[j];
			_At->rows_index[j][jc] = i;
			_At->rows_values[j][jc] = A->rows_values[i][q];
			jcount[j] = jc + 1;
		}
	}
	nb_soft_free_mem(memsize, jcount);
}

void nb_sparse_transpose(nb_sparse_t *A)
{
	for (uint32_t i = 0; i < A->N; i++) {
		for (uint32_t q = 0; A->rows_index[i][q] < i; q++) {
			uint32_t j = A->rows_index[i][q];
			uint32_t jc = nb_sparse_bsearch_row(A, j, i, 0,
							 A->rows_size[j]-1);
			double aux = A->rows_values[i][q];
			A->rows_values[i][q] = A->rows_values[j][jc];
			A->rows_values[j][jc] = aux;
		}
	}
}

uint32_t nb_sparse_get_size(const nb_sparse_t *const A)
{
	return A->N;
}


uint32_t nb_sparse_get_nnz(const nb_sparse_t *const A)
{
	uint32_t nnz = 0;
	for (uint32_t i = 0; i < A->N; i++)
		nnz += A->rows_size[i];
	return nnz;
}

double nb_sparse_get_asym(const nb_sparse_t *const A)
{
	double sum = 0;
	for (uint32_t i = 0; i < A->N; i++) {
		for (uint32_t j = 0; A->rows_index[i][j] < i; j++) {
			double Aij = A->rows_values[i][j];
			uint32_t k = A->rows_index[i][j];
			uint32_t l = nb_sparse_bsearch_row(A, k, i, 0,
							A->rows_size[k]-1);
			double Aji = 0.0;
			if (A->rows_index[k][l] == i)
				Aji = A->rows_values[k][l];

			sum += POW2(Aij - Aji); 
		}
	}
	return sqrt(sum);
}

double nb_sparse_get_frobenius_norm(const nb_sparse_t *const A)
{
	double sum = 0;
	for (uint32_t i = 0; i < A->N; i++) {
		uint32_t N_cols = A->rows_size[i];
		for (uint32_t j = 0; j < N_cols; j++) {
			double Aij = A->rows_values[i][j];
			sum += POW2(Aij); 
		}
	}
	return sqrt(sum);
}

void nb_sparse_multiply_scalar(nb_sparse_t* A, double scalar,
				uint32_t omp_parallel_threads)
{
#pragma omp parallel for num_threads(omp_parallel_threads) schedule(guided)
	for (uint32_t i = 0; i < A->N; i++) {
		for (uint32_t j = 0; j < A->rows_size[i]; j++)
			A->rows_values[i][j] *= scalar;
	}

}

void nb_sparse_multiply_vector(const nb_sparse_t* A, const double* in,
				double* out, uint32_t omp_parallel_threads)
{
#pragma omp parallel for num_threads(omp_parallel_threads) schedule(guided)
	for (uint32_t i = 0; i < A->N; i++) {
		out[i] = 0;
		for (uint32_t j = 0; j < A->rows_size[i]; j++)
			out[i] += A->rows_values[i][j] * in[A->rows_index[i][j]];
	}
}

void nb_sparse_set_Dirichlet_condition(nb_sparse_t* A, double* RHS,
					uint32_t idx, double value)
{
	for (uint32_t j = 0; j < A->rows_size[idx]; j++){
		uint32_t jdx = A->rows_index[idx][j];
		if (idx == jdx) {
			A->rows_values[idx][j] = 1.0;
			RHS[idx] = value;
		} else {
			A->rows_values[idx][j] = 0.0;
			double var = nb_sparse_get_and_set(A, jdx, idx, 0.0);
			RHS[jdx] -= var * value;
		}
	}
}
