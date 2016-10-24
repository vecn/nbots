#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "nb/math_bot.h"
#include "nb/memory_bot.h"
#include "nb/container_bot.h"
#include "nb/solver_bot/sparse/sparse.h"

#include "sparse_struct.h"

#define POW2(a) ((a)*(a))

static inline double get_norm(double* x, uint32_t N);

static int meta_compare_data_bycol(const void *a, const void *b);

static inline double get_norm(double* x, uint32_t N)
{
	double n = 0;
	for(uint32_t i=0; i<N; i++)
		n += x[i]*x[i];
	return sqrt(n);
}

nb_sparse_t* nb_sparse_create(const nb_graph_t *const restrict graph,
				const uint32_t *const restrict perm,
				uint32_t vars_per_node)
{
	nb_sparse_t* A = sparse_allocate(graph->N * vars_per_node);

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

int nb_sparse_solve_Gauss_Seidel
(const nb_sparse_t *const A, 
 const double *const b,
 double *_x,                /* Out */
 uint32_t max_iter, double tolerance,
 uint32_t* niter_performed,     /* Out (NULL if not required) */
 double* tolerance_reached, /* Out (NULL if not required) */
 uint32_t omp_parallel_threads){
	/* Allocate RHS vector */
	double* c = (double*)nb_allocate_mem(nb_sparse_get_size(A)*sizeof(double));
	/* Start iterations */
	register uint32_t k = 0;
	double error = 1e10;
	while(k < max_iter){
		/* Generate RHS vector */
		error = 0;
#pragma omp parallel for reduction(+:error) num_threads(omp_parallel_threads)
		for(uint32_t i=0; i < nb_sparse_get_size(A); i++){
			c[i] = b[i];
			double error_i = b[i];
			for(uint32_t j=0; j < A->rows_size[i]; j++){
				error_i -= A->rows_values[i][j] * _x[A->rows_index[i][j]];
				/* Use only the Upper triangular matrix to create the RHS */
				if(A->rows_index[i][j] > i)
					c[i] -= A->rows_values[i][j] * _x[A->rows_index[i][j]];
			}
			error += POW2(error_i);
		}
		/* Check convergence */
		if(error <= tolerance*tolerance)
			break;

		/* Solve by forward substitution */
		nb_sparse_forward_solve(A, c, _x);

		/* Increase iterations */
		k++;
	}
	/* Free memory */
	nb_free_mem(c);

	/* Successful exit */
	if(niter_performed != NULL) *niter_performed = k;

	if(tolerance_reached != NULL)
		*tolerance_reached = sqrt(error);
  
	if(error > tolerance*tolerance)
		return 1;
	return 0;
}

int nb_sparse_solve_conjugate_gradient
(const nb_sparse_t *const A, 
 const double *const b, 
 double *_x,                /* Out */
 uint32_t max_iter, double tolerance,
 uint32_t* niter_performed,     /* Out (NULL if not required) */
 double* tolerance_reached, /* Out (NULL if not required) */
 uint32_t omp_parallel_threads)
/* Return the num of iterations */
{
	/* Solve Ax = b with Conjugate Gradient method */
	double *g = nb_allocate_zero_mem(A->N * sizeof(double));
	double *p = nb_allocate_zero_mem(A->N * sizeof(double));
	double *w = nb_allocate_zero_mem(A->N * sizeof(double));
	double dot_gg = 0;

#pragma omp parallel for reduction(+:dot_gg) num_threads(omp_parallel_threads) schedule(guided)
	for(uint32_t i=0; i< A->N; i++){
		double sum = 0;
		for(uint32_t j=0; j< A->rows_size[i]; j++)
			sum += A->rows_values[i][j] * _x[A->rows_index[i][j]];
		g[i] = sum - b[i];
		p[i] = -g[i];
		dot_gg += g[i]*g[i];
	}
	uint32_t k = 0;
	while(dot_gg > tolerance*tolerance && k < max_iter){
		double dot_pw = 0;
		dot_gg = 0;
#pragma omp parallel for reduction(+:dot_pw, dot_gg) num_threads(omp_parallel_threads) schedule(guided)
		for(uint32_t i = 0; i< A->N; i++){
			w[i] = 0;
			for(uint32_t j = 0; j< A->rows_size[i]; j++)
				w[i] += A->rows_values[i][j] * p[A->rows_index[i][j]];
			dot_pw += p[i]*w[i];
			dot_gg += g[i]*g[i];
		}
		double alphak = dot_gg/dot_pw;
		double dot_gkgk = 0;
#pragma omp parallel for reduction(+:dot_gkgk) num_threads(omp_parallel_threads) schedule(guided)
		for(uint32_t i=0; i< A->N; i++){
			_x[i] += alphak*p[i];
			g[i] += alphak*w[i];
			dot_gkgk += g[i]*g[i];
		}
		double betak = dot_gkgk/dot_gg;
#pragma omp parallel for num_threads(omp_parallel_threads)
		for(uint32_t i=0; i< A->N; i++)
			p[i] = -g[i] + betak * p[i];
		k++;
	}
	/* Free memory */
	nb_free_mem(g);
	nb_free_mem(p);
	nb_free_mem(w);
  
	if(niter_performed != NULL) niter_performed[0]= k;

	if(tolerance_reached != NULL) *tolerance_reached = sqrt(dot_gg);

	if(dot_gg > tolerance*tolerance)
		return 1;

	return 0;
}

void nb_sparse_eigen_power(const nb_sparse_t* const A, int h,
			    double **_eigenvecs,/* Out */ 
			    double *_eigenvals, /* Out */
			    int *it,            /* Out */
			    double tolerance,
			    uint32_t omp_parallel_threads){
	/* The program must receive all the pointers allocated, where
	 *  > A is a nb_sparse_t matrix
	 *  > _eigenvecs is an array of size h to store h eigenvectors.
	 *  > _eigenvals is an array of size h to store the h greatest
	 *    eigenvalues approximated.
	 *  > h is the number of eigenvalues to be computed.
	 *  > '*it' will store (after computation) the iterations needed 
	 *    to compute each eigenvalue (is a return value).
	 */
		
	/* Declare structures and variables to be used */
	uint32_t i, j, c, d; /* Iterative variables */
	double pnorm, rnorm2;
	/* Allocate memory for structures */
	double *p = nb_allocate_zero_mem(A->N * sizeof(double));

	/* Deflation power method */
	for (i=0; i<h; i++){
		it[i] = 0;
		rnorm2 = 1;
		/* Initialize q0 such that ||qk||=1 */
		_eigenvecs[i][0] = 1;
		for(j=1; j < A->N; j++)
			_eigenvecs[i][j] = 0;

		for(c=0; c < A->N; c++){
			p[c] = 0;
			if(A->rows_index[c][0] == 0)
				p[c] = A->rows_values[c][0];
		}
		/* Start loop */
		while(rnorm2 > POW2(tolerance)){
			/* Step 1 */
			pnorm = get_norm(p, A->N);
			for(c=0; c < A->N; c++)
				_eigenvecs[i][c] = p[c]/pnorm;
			/* Step 2 */
			for(j=0; j<i; j++){
				double alpha = 0;

#pragma omp parallel for reduction(+:alpha) num_threads(omp_parallel_threads) schedule(guided) private(c)
				for(c=0; c < A->N; c++)
					alpha += _eigenvecs[i][c]*_eigenvecs[j][c];

#pragma omp parallel for num_threads(omp_parallel_threads) private(c)
				for(c=0; c < A->N; c++)
					_eigenvecs[i][c] -= alpha*_eigenvecs[j][c];
			}
			/* Step 3 */
			/* Paralelize the operation pk = A*qk */

#pragma omp parallel for schedule(guided) num_threads(omp_parallel_threads) private(c, d)
			for(c=0; c < A->N; c++){
				p[c] = 0;
				for(d=0; d < A->rows_size[c]; d++)
					p[c] += A->rows_values[c][d]
						*_eigenvecs[i][A->rows_index[c][d]];
			}
			/* Step 4 */
			double lambda = 0;
			for(c=0; c < A->N; c++)
				lambda += _eigenvecs[i][c]*p[c];
			_eigenvals[i] = lambda;
			/* Step 5 and 6 */
			rnorm2 = 0;
			for(c=0; c < A->N; c++)
				rnorm2 += POW2(p[c]-lambda*_eigenvecs[i][c]);
			it[i]++;
		}
	}

	/* Free memory */
	nb_free_mem(p);
}

int nb_sparse_eigen_ipower(const nb_sparse_t *const A,
			    nb_solver_t solver,
			    int h, double mu, 
			    double **_eigenvecs,/* Out */ 
			    double *_eigenvals, /* Out */
			    int* it,            /* Out */
			    double tolerance,
			    uint32_t omp_parallel_threads)
{
	/* The program must receive all the pointers allocated, where
	 *  > A is a nb_sparse_t matrix
	 *  > _eigenvecs is an array of size h to store h eigenvectors.
	 *  > _eigenvals is an array of size h to store the h
	 *    eigenvalues approximated.
	 *  > h is the number of eigenvalues to be computed.
	 *  > '*it' will store (after computation) the iterations needed 
	 *    to compute each eigenvalue (is a return value).
	 */

	/* Declare structures and variables to be used */
	uint32_t i, j, c, d; /* Iterative variables */
	double pnorm, rnorm2;
	/* Allocate memory for structures */
	double* p = nb_allocate_zero_mem(A->N * sizeof(double));
	double* z = nb_allocate_zero_mem(A->N * sizeof(double));

	/* Set M in A (copy ptr to modify const A, it will be restored) */
	nb_sparse_t* A_ptr_copy = (nb_sparse_t*)A;
	if (mu != 0.0)
		for (c = 0; c < A->N; c++)
			nb_sparse_add(A_ptr_copy, c, c, -mu);            /* M = A - mu*I */

	/* LU Decomposition in case of LU Solver */
	nb_sparse_t *L = NULL;
	nb_sparse_t *U = NULL;
	if (NB_SOLVER_CHK == solver) {
		nb_sparse_alloc_LU(A, &L, &U);
		nb_sparse_decompose_Cholesky(A, L, U, omp_parallel_threads);
	} else if (NB_SOLVER_LUD == solver) {
		nb_sparse_alloc_LU(A, &L, &U);
		nb_sparse_decompose_LU(A, L, U, omp_parallel_threads);
	}

	/* Deflation inverse power method */
	for (i = 0; i < h; i++) {
		it[i] = 0;
		rnorm2 = 1;
		/* Initialize q0 such that ||qk||=1 */
		_eigenvecs[i][0] = 1;
		for(j=1; j< A->N; j++)
			_eigenvecs[i][j] = 0;

		/* Start loop */
		double rnorm2_diff = 1;
		while(rnorm2 > POW2(tolerance) && rnorm2_diff > POW2(tolerance)){
			/* Step 1 */
			if (NB_SOLVER_CHK == solver ||
			    NB_SOLVER_LUD == solver)
				nb_sparse_solve_LU(L, U, _eigenvecs[i], p);
			else if (NB_SOLVER_CGJ == solver)
				nb_sparse_solve_CG_precond_Jacobi(A,_eigenvecs[i], p, 
								   nb_sparse_get_size(A)*10, 1e-3, 
								   NULL, NULL,
								   omp_parallel_threads);
			else
				nb_sparse_solve_conjugate_gradient(A,_eigenvecs[i], p, 
								    nb_sparse_get_size(A)*10, 1e-3, 
								    NULL, NULL,
								    omp_parallel_threads);
			/* Step 2 */
			pnorm = get_norm(p, A->N);
			for(c=0; c < A->N; c++)
				_eigenvecs[i][c] = p[c]/pnorm;
			/* Step 3 */
			for (j = 0; j < i; j++){
				double alpha = 0;

#pragma omp parallel for reduction(+:alpha) num_threads(omp_parallel_threads) schedule(guided) private(c)
				for(c=0; c < A->N; c++)
					alpha += _eigenvecs[i][c]*_eigenvecs[j][c];

#pragma omp parallel for num_threads(omp_parallel_threads) private(c)
				for(c=0; c < A->N; c++)
					_eigenvecs[i][c] -= alpha*_eigenvecs[j][c];
			}
			/* Step 4 */
			/* Paralelize the operation zk = A*qk */
#pragma omp parallel for num_threads(omp_parallel_threads) schedule(guided) private(d) private(c)
			for(c=0; c < A->N; c++){
				z[c] = 0;
				for(d=0; d < A->rows_size[c]; d++){
					double aii = A->rows_values[c][d];
					if(c == d) aii += mu;
					z[c] += aii
						*(_eigenvecs[i][A->rows_index[c][d]]);
				}
			}
			/* Step 5 */
			double sigma = 0;
			for(c=0; c < A->N; c++)
				sigma += _eigenvecs[i][c]*z[c];
			_eigenvals[i] = sigma;
			/* Step 6 and 7 */
			rnorm2_diff = rnorm2;
			rnorm2 = 0;
			for(c=0; c < A->N; c++)
				rnorm2 += POW2(z[c]-sigma*_eigenvecs[i][c]);
			rnorm2_diff = fabs(rnorm2_diff-rnorm2);
			it[i]++;
		}
	}
	/* Restore A */
	if (mu != 0.0)
		for (c = 0; c < A->N; c++)
			nb_sparse_add(A_ptr_copy, c, c, mu);  /* A = M + mu*I */
    
	/* Destroy LU decomposition */
	if (NB_SOLVER_CHK == solver || NB_SOLVER_LUD == solver) {
		nb_sparse_destroy(U);
		nb_sparse_destroy(L);
	}
	/* Free memory */
	nb_free_mem(p);
	nb_free_mem(z);

	return 0;
}


void nb_sparse_eigen_lanczos(const nb_sparse_t* const A,
			      double *_eigenmax,/* Out */ 
			      double *_eigenmin,/* Out */
			      int* it,          /* Out */
			      double tolerance,
			      uint32_t omp_parallel_threads)
{
	/* The program must receive all the pointers allocated, where
	 *  > A is a nb_sparse_t matrix
	 *  > _eigenvecs is an array of size h to store h eigenvectors.
	 *  > _eigenvals is an array of size h to store the h greatest
	 *    eigenvalues approximated.
	 *  > h is the number of eigenvalues to be computed.
	 *  > '*it' will store (after computation) the iterations needed 
	 *    to compute each eigenvalue (is a return value).
	 */
	/* Declare structures and variables to be used */
	int i, j;

	/* Declare structures and variables to be used */
	double* alpha = nb_allocate_zero_mem(A->N * sizeof(double));
	double* beta = nb_allocate_zero_mem(A->N * sizeof(double));
	double* v = nb_allocate_zero_mem(A->N * sizeof(double));
	double* w = nb_allocate_zero_mem(A->N * sizeof(double));
	double* v_prev = nb_allocate_zero_mem(A->N * sizeof(double));

	*it = 0;
	beta[*it] = 0;
	v[0] = 1;
	/* Start iterations */
	double normw2 = 1;
	double delta_eigen = 1;

	while (delta_eigen >= POW2(tolerance) &&
	       normw2 >= tolerance && *it < A->N) {
		/* Step 1 and 2 */
		double ak = 0;
#pragma omp parallel for reduction(+:ak) num_threads(omp_parallel_threads) private(i, j)
		for (i = 0; i < A->N; i++) {
			w[i] = 0;
			for (j = 0; j < A->rows_size[i]; j++)
				w[i] += A->rows_values[i][j]*v[A->rows_index[i][j]];
			w[i] -= beta[*it]*v_prev[i];
			ak += w[i]*v[i];
		}
		alpha[*it] = ak;
		/* Step 3 and 4 */
		*it = *it+1;
		normw2 = 0;
		for (i = 0; i < A->N; i++) {
			w[i] -= ak*v[i];
			normw2 += w[i]*w[i];
		}
		normw2 = sqrt(normw2);
		beta[*it] = normw2;
		/* Step 5 */
		for (i = 0; i < A->N; i++) {
			v_prev[i] = v[i];
			v[i] = w[i]/beta[*it];
		}
		/* Step 6 and 7 */
		if (*it > 1) {
			double delta_eig1 = *_eigenmax;
			double delta_eigk = *_eigenmin;
			nb_sparse_eigen_givens(alpha, beta, 1, _eigenmax, tolerance, *it);
			nb_sparse_eigen_givens(alpha, beta, *it, _eigenmin, tolerance, *it);
      
			delta_eigen = fabs(delta_eig1-*_eigenmax);
			double tmp = fabs(delta_eigk-*_eigenmin);
			if (tmp > delta_eigen)
				delta_eigen = tmp;
		} else {
			*_eigenmax = ak;
			*_eigenmin = ak;
		}
	}

	/* Free memory */
	nb_free_mem(alpha);
	nb_free_mem(beta);
	nb_free_mem(w);
	nb_free_mem(v);
	nb_free_mem(v_prev);
}

void nb_sparse_eigen_givens(const double* const main_diag, 
			     const double* const uplw_diag,
			     int i, double *_eigenvalue,
			     double tolerance, uint32_t N)
{
	/* Compute the eigen values of a symmetric tridiagonal 
	 * matrix (using Givens method), where
	 *   > *main_diag is the main diagonal of the matrix.
	 *   > *uplw_diag is the upper and lower diagonal, the
	 *     first value will be ignored.
	 *   > 'i' is the position of the eigen value to compute, 
	 *     assuming:
	 *        l1 > l2 > l2 > ... > ln
	 *     where li is the ith eigenvalue.
	 *   > The program return through '*_eigenvalue' the
	 *     resulting value.
	 */
	/* Iterative variables */
	int l;
	/* Intitialize and allocate structures */
	double* p = nb_allocate_zero_mem((N+1) * sizeof(double));

	/* Init algorithm */
	int k = N;
	double a = main_diag[0]-fabs(uplw_diag[1]);
	double tmp = main_diag[k-1]-fabs(uplw_diag[k-1]);
	if (tmp < a)
		a = tmp;
	double b = main_diag[0]+fabs(uplw_diag[1]);
	tmp = main_diag[k-1]+fabs(uplw_diag[k-1]);
	if (tmp > b)
		b = tmp;
	for (l=1; l<k-1; l++) {
		tmp = main_diag[l]-
			fabs(uplw_diag[l+1])-fabs(uplw_diag[l]);
		if (tmp < a)
			a = tmp;
		tmp = main_diag[l]+
			fabs(uplw_diag[l+1])+fabs(uplw_diag[l]);
		if (tmp > b)
			b = tmp;
	}
	/* Init iterations */
	while (fabs(b-a) > (fabs(a)+fabs(b)) * tolerance) {
		/* Step 1 */
		*_eigenvalue = (a+b)/2;
		/* Step 2 */
		double r = 0;
		double s = 0;
		p[0] = 1;
		p[1] = main_diag[0] - *_eigenvalue;
		for (l = 1; l < k; l++) {
			p[l+1] = (main_diag[l]-*_eigenvalue)*p[l]-
				(uplw_diag[l])*(uplw_diag[l])*p[l-1];
		}
		for (l=1; l<k+1; l++) {
			if (p[l]*p[l-1] <= 0)
				r++;
			if (p[l] == 0)
				s++;
		}
		double gamma = r-s;
		/* Step 3 */
		if (gamma > k-i)
			b = *_eigenvalue;
		else
			a = *_eigenvalue;
	}

	/* Free memory */
	nb_free_mem(p);
}

static int meta_compare_data_bycol(const void *a, const void *b)
{
	if((*(double**)a)[0] == (*(double**)b)[0])
		return (*(double**)a)[1] - (*(double**)b)[1];
	else
		return (*(double**)a)[0] - (*(double**)b)[0];
}
