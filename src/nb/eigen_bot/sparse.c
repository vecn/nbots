/******************************************************************************
 *   Sparse Bot: Linear Algebra for sparse and symmetric matrices.            *
 *   2011-2015 Victor Eduardo Cardoso Nungaray                                *
 *   Twitter: @victore_cardoso                                                *
 *   email: victorc@cimat.mx                                                  *
 ******************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include "nb/math_bot.h"
#include "nb/container_bot/array.h"
#include "nb/eigen_bot/sparse.h"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define POW2(a) ((a)*(a))

static inline double get_norm(double* x, uint32_t N);

static int meta_compare_data_bycol(const void *a, const void *b);

static inline vcn_sparse_t* sparse_allocate(uint32_t N);

static void cholesky_symbolic(const vcn_sparse_t *const A,
			      vcn_sparse_t *_L, vcn_sparse_t* _Lt, 
			      uint32_t ktrunc);

static inline uint32_t sparse_bsearch_row(const vcn_sparse_t *const A,
					  uint32_t i, uint32_t col, 
					  int imin, int imax);
static inline void matrix_r_solver(const double *const A, int n, 
				   double *d, double *b);

struct vcn_sparse_s {
	double **rows_values;
	uint32_t **rows_index;
	uint32_t *rows_size;
	uint32_t N;
} ;

static inline double get_norm(double* x, uint32_t N)
{
	double n = 0;
	for(uint32_t i=0; i<N; i++)
		n += x[i]*x[i];
	return sqrt(n);
}

static int meta_compare_data_bycol(const void *a, const void *b)
{
	if((*(double**)a)[0] == (*(double**)b)[0])
		return (*(double**)a)[1] - (*(double**)b)[1];
	else
		return (*(double**)a)[0] - (*(double**)b)[0];
}

/* Warning: The methods of sparse struct doesn't have any 
 *  validation, because they are thought to be fast.
 *  Ensure all pointers are allocated and initialized, and
 *  all index are inside bounds, etc.
 */

static inline vcn_sparse_t* sparse_allocate(uint32_t N)
{
	vcn_sparse_t* A = malloc(sizeof(*A));
	A->rows_values = malloc(N * sizeof(*(A->rows_values)));
	A->rows_index = malloc(N * sizeof(*(A->rows_index)));
	A->rows_size = calloc(N, sizeof(*(A->rows_size)));
	A->N = N;
	return A;
}

vcn_sparse_t* vcn_sparse_create(const vcn_graph_t *const restrict graph,
				const uint32_t *const restrict perm,
				uint32_t vars_per_node)
{
	vcn_sparse_t* A = sparse_allocate(graph->N * vars_per_node);

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
				calloc(row_size,
				       sizeof(*(A->rows_values[irow])));
			A->rows_index[irow] = 
				malloc(row_size *
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
			vcn_qsort(A->rows_index[irow], A->rows_size[irow], 
				  sizeof(uint32_t), vcn_compare_uint32);
		}
	}
	return A;
}

vcn_sparse_t* vcn_sparse_clone(vcn_sparse_t* A)
{
	vcn_sparse_t* Acopy = sparse_allocate(A->N);
	for (uint32_t i = 0; i < A->N; i++) {
		Acopy->rows_index[i] = calloc(A->rows_size[i], 
					      sizeof(**(Acopy->rows_index)));
		memcpy(Acopy->rows_index[i], A->rows_index[i], 
		       A->rows_size[i] * sizeof(**(Acopy->rows_index)));

		Acopy->rows_values[i] = calloc(A->rows_size[i],
					       sizeof(**(Acopy->rows_values)));
		memcpy(Acopy->rows_values[i], A->rows_values[i],
		       A->rows_size[i] * sizeof(**(Acopy->rows_values)));

		Acopy->rows_size[i] = A->rows_size[i];
	}
	return Acopy;
}

void vcn_sparse_save(const vcn_sparse_t *const A, const char* filename)
{
	FILE *fp = fopen(filename, "w");
	if (fp == NULL) {
		printf("Error: Opening the file %s to save the matrix.\n", filename);
		return;
	}
	fprintf(fp, "%i %i\n", A->N, A->N);
	for (uint32_t i=0; i < A->N; i++) {
		for (uint32_t j=0; j < A->rows_size[i]; j++)
			fprintf(fp, "%i %i %e \n", i, A->rows_index[i][j], A->rows_values[i][j]);
	}
	fclose(fp);
}

void vcn_sparse_destroy(vcn_sparse_t* A)
{
	/* Clear all rows */
	for (uint32_t i = 0; i < A->N; i++) {
		free(A->rows_values[i]);
		free(A->rows_index[i]);
	}
	free(A->rows_values);
	free(A->rows_index);
	free(A->rows_size);
	free(A);
}

static inline uint32_t sparse_bsearch_row(const vcn_sparse_t *const A,
					  uint32_t i, uint32_t col, 
					  int imin, int imax)
{
	uint32_t index = (imin+imax)/2;
	while (A->rows_index[i][index] != col && imin <= imax) {
		index = (imin+imax)/2;
		if (A->rows_index[i][index] < col)
			imin = index+1;
		else
			imax = index-1;
	}
	return index;
}

void vcn_sparse_reset(vcn_sparse_t *A)
{
	for(uint32_t i = 0; i < A->N; i++)
		memset(A->rows_values[i], 0,
		       A->rows_size[i] * sizeof(*(A->rows_values)));  
}

void vcn_sparse_set(vcn_sparse_t *A, uint32_t i, uint32_t j, double value)
{
	uint32_t index = sparse_bsearch_row(A, i, j, 0, A->rows_size[i]-1);
	if (A->rows_index[i][index] == j) {
		A->rows_values[i][index] = value;
	} else {
		/* OPPORTUNITY: Send this to some logfile */
		printf("ERROR: Entry of sparse matrix is not allocated.\n");
		printf("    -> vcn_sparse_set(*A=%p, i=%i, j=%i, value=%lf)\n",
		       (void*) A, i, j, value);
		exit(1);
	}
}

void vcn_sparse_set_identity_row(vcn_sparse_t* A, uint32_t row)
{
	memset(A->rows_values[row], 0,
	       A->rows_size[row] * sizeof(**(A->rows_values)));
	vcn_sparse_set(A, row, row, 1.0);
}

void vcn_sparse_make_diagonal(vcn_sparse_t* A, double diag_val)
{
	for (uint32_t i = 0; i < A->N; i++){
		memset(A->rows_values[i], 0, 
		       A->rows_size[i] * sizeof(**(A->rows_values)));
		vcn_sparse_set(A, i, i, diag_val);
	}
}

double vcn_sparse_get(const vcn_sparse_t *const A, uint32_t i, uint32_t j)
{
	uint32_t index = sparse_bsearch_row(A, i, j, 0, A->rows_size[i]-1);
	if (A->rows_index[i][index] == j)
		return A->rows_values[i][index];
	else
		return 0.0;
}

double vcn_sparse_get_and_set(vcn_sparse_t *A, uint32_t i, uint32_t j, double value)
{
	uint32_t index = sparse_bsearch_row(A, i, j, 0, A->rows_size[i]-1);
	if (A->rows_index[i][index] == j) {
		double val = A->rows_values[i][index];
		A->rows_values[i][index] = value;
		return val;
	} else {
		/* OPPORTUNITY: Send this to some logfile */
		printf("ERROR: Entry of sparse matrix is not allocated.\n");
		printf("    -> vcn_sparse_get_and_set(*A=%p, i=%i, j=%i, value=%lf)\n",
		       (void*) A, i, j, value);
		exit(1);
	}
}

bool vcn_sparse_is_non_zero(const vcn_sparse_t *const A, uint32_t i, uint32_t j)
{
	uint32_t index = sparse_bsearch_row(A, i, j, 0, A->rows_size[i]-1);
	if(A->rows_index[i][index] == j)
		return true;
	else
		return false;
}

uint32_t vcn_sparse_memory_used(const vcn_sparse_t *const A)
{
	uint32_t size = sizeof(short) + 3 * sizeof(uint32_t) + 3*sizeof(void*);
	size += A->N * (2*sizeof(void*) + sizeof(uint32_t));
	for (uint32_t i = 0; i < A->N; i++)
		size += A->rows_size[i] * (sizeof(uint32_t) + sizeof(double));
	return size;
}

void vcn_sparse_add(vcn_sparse_t *A, uint32_t i, uint32_t j, double value)
{
	uint32_t index = sparse_bsearch_row(A, i, j, 0, A->rows_size[i]-1);
	if (A->rows_index[i][index] == j) {
		A->rows_values[i][index] += value;
	} else {
		/* OPPORTUNITY: Send this to some logfile */
		printf("ERROR: Entry of sparse matrix is not allocated.\n");
		printf("    -> vcn_sparse_add(*A=%p, i=%i, j=%i, value=%lf)\n",
		       (void*) A, i, j, value);
		exit(1);
	}
}

void vcn_sparse_scale(vcn_sparse_t *A, double factor)
{
	for (uint32_t i = 0; i < A->N; i++) {
		for (uint32_t j = 0; j < A->rows_size[i]; j++)
			A->rows_values[i][j] *= factor;
	}
}
 
vcn_sparse_t* vcn_sparse_create_permutation
(const vcn_sparse_t *const A,
 const uint32_t *const perm,
 const uint32_t *const iperm)
{
	/* Use the vectors perm and iperm to compute Ar = PAP'. */ 
	vcn_sparse_t* _Ar = sparse_allocate(A->N);
	for (uint32_t i = 0; i < A->N; i++) {
		uint32_t j = perm[i];
		_Ar->rows_size[i] = A->rows_size[j];
		_Ar->rows_index[i] = calloc(_Ar->rows_size[i],sizeof(**(_Ar->rows_index)));
		_Ar->rows_values[i] = calloc(_Ar->rows_size[i],sizeof(**(_Ar->rows_values)));
		double** data2sort = malloc(_Ar->rows_size[i]*sizeof(*data2sort));
		for (uint32_t k = 0; k < A->rows_size[j]; k++) {
			uint32_t m = iperm[A->rows_index[j][k]];
			data2sort[k] = calloc(2, sizeof(**data2sort));
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
			free(data2sort[k]);
		free(data2sort);
	}
	return _Ar;
}

void vcn_sparse_fill_permutation(const vcn_sparse_t *const A, 
				 vcn_sparse_t* Ar,
				 const uint32_t *const perm,
				 const uint32_t *const iperm)
{
	/* Use the vectors perm and iperm to compute Ar = PAP'. */
	for (uint32_t i=0; i < A->N; i++) {
		uint32_t j = perm[i];
		for (uint32_t k = 0; k < A->rows_size[j]; k++) {
			uint32_t m = iperm[A->rows_index[j][k]];
			uint32_t idx = sparse_bsearch_row(Ar, i, m, 0, Ar->rows_size[i]-1);
			Ar->rows_values[i][idx] = A->rows_values[j][k];
		}
	}
}

double* vcn_sparse_create_vector_permutation
(const double *const b, 
 const uint32_t *const perm,
 uint32_t N){
	double* br = (double*)malloc(N * sizeof(double));
	for(uint32_t i=0; i < N; i++)
		br[i] = b[perm[i]];
	return br;
}

void vcn_sparse_transpose(vcn_sparse_t *A, vcn_sparse_t *_At)
/* _At must be allocated and initialized */
{
	uint32_t *jcount = calloc(A->N,sizeof(*jcount));
	for (int i = 0; i < A->N; i++) {
		for (int q = 0; q < A->rows_size[i]; q++) {
			uint32_t j = A->rows_index[i][q];
			_At->rows_index[j][jcount[j]] = i;
			_At->rows_values[j][jcount[j]++] = A->rows_values[i][q];
		}
	}
	free(jcount);
}

uint32_t vcn_sparse_get_size(const vcn_sparse_t *const A)
{
	return A->N;
}


uint32_t vcn_sparse_get_nnz(const vcn_sparse_t *const A)
{
	uint32_t nnz = 0;
	for (uint32_t i = 0; i < A->N; i++)
		nnz += A->rows_size[i];
	return nnz;
}

void vcn_sparse_multiply_scalar(vcn_sparse_t* A, double scalar,
				uint32_t omp_parallel_threads)
{
#pragma omp parallel for num_threads(omp_parallel_threads) schedule(guided)
	for (uint32_t i = 0; i < A->N; i++) {
		for (uint32_t j = 0; j < A->rows_size[i]; j++)
			A->rows_values[i][j] *= scalar;
	}

}

void vcn_sparse_multiply_vector(vcn_sparse_t* A, double* in, double* out,
				uint32_t omp_parallel_threads)
{
#pragma omp parallel for num_threads(omp_parallel_threads) schedule(guided)
	for (uint32_t i = 0; i < A->N; i++) {
		out[i] = 0;
		for (uint32_t j = 0; j < A->rows_size[i]; j++)
			out[i] += A->rows_values[i][j] * in[A->rows_index[i][j]];
	}
}

int vcn_sparse_spy_plot_as_png(const vcn_sparse_t *const A,
			       const char* url, uint32_t size,
			       bool enable_zeros_allocated,
			       bool enable_color)
{
//	/* Spyplot properties */
//	int border = 3;
//	int rgb_border[] = {120, 120, 120};
//	int s_palete = 4;
//	float p_palete[] = {0, 0.0001, 0.05, 0.2};
//	int rgb_palete[12];
//	if(enable_color){
//		int rgb [] =
//			{0, 0, 0,        /* Black */
//			 99, 0, 220,     /* Purple */
//			 255, 60, 0,     /* Orange */
//			 255, 255, 0};   /* Yellow */
//		memcpy(rgb_palete, rgb, 12 * sizeof(int));
//	}else{
//		int rgb[] = 
//			{255, 255, 255,    /* White */
//			 200, 200, 200,    /* Dark gray */
//			 50, 50, 50,       /* Gray */
//			 0, 0, 0};         /* Black */
//		memcpy(rgb_palete, rgb, 12 * sizeof(int));
//	}
//
//	/* Create spyplot */
//	int i, j, u, v;                   /* Iterative variables */
//	int plot_size = size-2*border;
//	float cell_size = ((float)A->N)/plot_size;
//	float cell_max = 0;
//	float* cells = (float*)calloc(plot_size*plot_size, sizeof(float));
//	for(i=0; i<plot_size; i++){
//		for(j=0; j<plot_size; j++){
//			/* Start iterations over matrix */
//			for(u=floor(i*cell_size); u<ceil((i+1)*cell_size) && u < A->N; u++){
//				v = 0;
//				while(v < A->rows_size[u] && 
//				      A->rows_index[u][v] < ceil((j+1)*cell_size)){
//					if(!enable_zeros_allocated){
//						if(A->rows_values[u][v] == 0.0){
//							v++;
//							continue;
//						}
//					}
//					if(A->rows_index[u][v] >= floor(j*cell_size)){
//						float value = 1;
//						if(A->N > plot_size){
//							if(u == floor(i*cell_size))
//								value *= ceil(i*cell_size)-i*cell_size;
//							else if(u == ceil((i+1)*cell_size))
//								value *= (i+1)*cell_size-floor((i+1)*cell_size);
//	    
//							if(A->rows_index[u][v] == floor(j*cell_size))
//								value *= ceil(j*cell_size)-j*cell_size;
//							else if(A->rows_index[u][v] == ceil((j+1)*cell_size)-1)
//								value *= (j+1)*cell_size-floor((j+1)*cell_size);
//						}
//						cells[i*plot_size+j] += value;	      
//					}
//					v++;
//				}
//			}
//			/* Finish iterations over matrix */
//			if(cells[i*plot_size+j] > cell_max)
//				cell_max = cells[i*plot_size+j];
//		}      
//	}
//
//	/* Normalize cells matrix */
//	for(i=0; i<plot_size*plot_size; i++)
//		cells[i] /= cell_max;
//
//	/* Declare structures and variables to be used */
//	png_structp png_ptr = NULL;
//	png_infop info_ptr = NULL;
//	png_byte** row_pointers = NULL;
//
//	/* Image Properties */
//	int pixel_size = 3;
//	int depth = 8;
//	uint w = size;
//	uint h = size;
//  
//	/* Open url to write PNG file */
//	FILE *fp = fopen(url, "wb");
//
//	/* Verify all pointers are working */
//	if(fp == NULL){
//		printf("ERROR: Open url to save PNG failed\n");
//		return 1;
//	}
//  
//	png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
//	if(png_ptr == NULL){
//		printf("ERROR: PNG create write struct failed.\n");
//		fclose(fp);
//		return 2;
//	}
//  
//	info_ptr = png_create_info_struct(png_ptr);
//	if(info_ptr == NULL){
//		printf("ERROR: PNG create info struct failed.\n");
//		png_destroy_write_struct(&png_ptr, &info_ptr);
//		fclose(fp);
//		return 3;
//	}
//  
//	if(setjmp(png_jmpbuf(png_ptr))){
//		printf("ERROR: PNG failure.\n");
//		png_destroy_write_struct(&png_ptr, &info_ptr);
//		fclose(fp);
//		return 4;
//	}
//
//	/* Set image attributes */
//	png_set_IHDR(png_ptr, info_ptr, w, h, depth, 
//		     PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
//		     PNG_COMPRESSION_TYPE_DEFAULT,
//		     PNG_FILTER_TYPE_DEFAULT);
//
//	/* Initialize pointers */
//	uint8_t r, g, b;
//	row_pointers = png_malloc(png_ptr, h*sizeof(png_byte*));
//	for(i=0; i<h; i++){
//		png_byte *row=
//			png_malloc(png_ptr, sizeof(uint8_t)*w*pixel_size);
//		row_pointers[i] = row;
//		for(j=0; j<w; j++){
//			/* Create bitmap */
//			if(j < border || j > size-border-1 ||
//			   i < border || i > size-border-1){
//				/* Create matrix borders */
//				r = rgb_border[0];
//				g = rgb_border[1];
//				b = rgb_border[2];
//			}else{
//				int x = j-border;
//				int y = i-border;
//				float val = cells[y*plot_size+x];
//				int p = 0;
//				while(p < s_palete && val > p_palete[p])
//					p++;
//				if(p == 0){
//					r = rgb_palete[0];
//					g = rgb_palete[1];
//					b = rgb_palete[2];
//				}else if(p == s_palete){
//					r = rgb_palete[(p-1)*3];
//					g = rgb_palete[(p-1)*3+1];
//					b = rgb_palete[(p-1)*3+2];
//				}else{
//					float val_scaled = (val-p_palete[p-1])/(p_palete[p]-p_palete[p-1]);
//					r = (rgb_palete[p*3]-rgb_palete[(p-1)*3])*val_scaled 
//						+ rgb_palete[(p-1)*3];
//					g = (rgb_palete[p*3+1]-rgb_palete[(p-1)*3+1])*val_scaled
//						+ rgb_palete[(p-1)*3+1];
//					b = (rgb_palete[p*3+2]-rgb_palete[(p-1)*3+2])*val_scaled
//						+ rgb_palete[(p-1)*3+2];
//				}
//			}
//			/* End bitmap creation */
//			*row++ = r;
//			*row++ = g;
//			*row++ = b;
//		}
//	}
//
//	/* Write the image data */
//	png_init_io(png_ptr, fp);
//	png_set_rows(png_ptr, info_ptr, row_pointers);
//	png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);
//
//	/* Free memory */
//	for(i=0; i<h; i++)
//		png_free(png_ptr, row_pointers[i]);
//	png_free(png_ptr, row_pointers);
//	png_destroy_write_struct(&png_ptr, &info_ptr);
//
//	free(cells);
//
//	/* Close url */
//	fclose(fp);
	return 0;
}

void vcn_sparse_set_Dirichlet_condition(vcn_sparse_t* A, double* RHS,
					uint32_t idx, double value)
{
	for (uint32_t j = 0; j < A->rows_size[idx]; j++){
		uint32_t jdx = A->rows_index[idx][j];
		if (idx == jdx) {
			A->rows_values[idx][j] = 1.0;
			RHS[idx] = value;
		} else {
			A->rows_values[idx][j] = 0.0;
			double var = vcn_sparse_get_and_set(A, jdx, idx, 0.0);
			RHS[jdx] -= var * value;
		}
	}
}

int vcn_sparse_solve_Gauss_Seidel
(const vcn_sparse_t *const A, 
 const double *const b,
 double *_x,                /* Out */
 uint32_t max_iter, double tolerance,
 uint32_t* niter_performed,     /* Out (NULL if not required) */
 double* tolerance_reached, /* Out (NULL if not required) */
 uint32_t omp_parallel_threads){
	/* Allocate RHS vector */
	double* c = (double*)malloc(vcn_sparse_get_size(A)*sizeof(double));
	/* Start iterations */
	register uint32_t k = 0;
	double error = 1e10;
	while(k < max_iter){
		/* Generate RHS vector */
		error = 0;
#pragma omp parallel for reduction(+:error) num_threads(omp_parallel_threads)
		for(uint32_t i=0; i < vcn_sparse_get_size(A); i++){
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
		vcn_sparse_forward_solve(A, c, _x);

		/* Increase iterations */
		k++;
	}
	/* Free memory */
	free(c);

	/* Successful exit */
	if(niter_performed != NULL) *niter_performed = k;

	if(tolerance_reached != NULL)
		*tolerance_reached = sqrt(error);
  
	if(error > tolerance*tolerance)
		return 1;
	return 0;
}

int vcn_sparse_solve_conjugate_gradient
(const vcn_sparse_t *const A, 
 const double *const b, 
 double *_x,                /* Out */
 uint32_t max_iter, double tolerance,
 uint32_t* niter_performed,     /* Out (NULL if not required) */
 double* tolerance_reached, /* Out (NULL if not required) */
 uint32_t omp_parallel_threads)
/* Return the num of iterations */
{
	/* Solve Ax = b with Conjugate Gradient method */
	double *g = (double*)calloc(A->N, sizeof(double));
	double *p = (double*)calloc(A->N, sizeof(double));
	double *w = (double*)calloc(A->N, sizeof(double));
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
	free(g);
	free(p);
	free(w);
  
	if(niter_performed != NULL) niter_performed[0]= k;

	if(tolerance_reached != NULL) *tolerance_reached = sqrt(dot_gg);

	if(dot_gg > tolerance*tolerance)
		return 1;

	return 0;
}

int vcn_sparse_solve_CG_precond_Jacobi
(const vcn_sparse_t *const A, 
 const double *const b,
 double *_x,                /* Out */
 uint32_t max_iter, double tolerance,
 uint32_t* niter_performed,     /* Out (NULL if not required) */
 double* tolerance_reached, /* Out (NULL if not required) */
 uint32_t omp_parallel_threads)
/* Return the num of iterations */
{
	/* Solve Ax = b with Conjugate Gradient preconditioned with jacobi */
	uint32_t vec_size = A->N * sizeof(double);
	char *memblock = calloc(1, 5 * vec_size);	
	double* g = (void*) memblock;
	double* p = (void*) (memblock + vec_size);
	double* q = (void*) (memblock + 2 * vec_size);
	double* w = (void*) (memblock + 3 * vec_size);
	double* Aii = (void*) (memblock + 4 * vec_size);
	double dot_gg = 0;

#pragma omp parallel for reduction(+:dot_gg) num_threads(omp_parallel_threads) schedule(guided)
	for (uint32_t i = 0; i < A->N; i++) {
		double sum = 0;
		for(uint32_t j=0; j< A->rows_size[i]; j++)
			sum += A->rows_values[i][j] * _x[A->rows_index[i][j]];
		g[i] = sum - b[i];
		Aii[i] = vcn_sparse_get(A,i,i);
		q[i] = g[i]/Aii[i];
		p[i] = -q[i];
		dot_gg += g[i]*g[i];
	}
	uint32_t k = 0;
	while (dot_gg > tolerance * tolerance && k < max_iter) {
		double dot_pw = 0;
		double dot_gq = 0;
		dot_gg = 0;

#pragma omp parallel for reduction(+:dot_pw, dot_gg, dot_gq) num_threads(omp_parallel_threads)
		for (uint32_t i = 0; i < A->N; i++){
			w[i] = 0;
			for (uint32_t j = 0; j <  A->rows_size[i]; j++)
				w[i] += A->rows_values[i][j] * p[A->rows_index[i][j]];
			dot_pw += p[i]*w[i];
			dot_gg += g[i]*g[i];
			dot_gq += g[i]*q[i];
		}
		double alphak = dot_gq/dot_pw;
		double dot_gkqk = 0;
		
#pragma omp parallel for reduction(+:dot_gkqk) num_threads(omp_parallel_threads)
		for(uint32_t i=0; i< A->N; i++){
			_x[i] += alphak*p[i];
			g[i] += alphak*w[i];
			q[i] = g[i]/Aii[i];
			dot_gkqk += g[i]*q[i];
		}

		double betak = dot_gkqk/dot_gq;
		
#pragma omp parallel for num_threads(omp_parallel_threads)
		for(uint32_t i=0; i< A->N; i++)
			p[i] = -q[i]+betak*p[i];
		k++;
	}
	/* Free memory */
	free(memblock);

	if (NULL != niter_performed)
		niter_performed[0]= k;

	if (NULL != tolerance_reached)
		*tolerance_reached = sqrt(dot_gg);

	if (dot_gg > tolerance*tolerance)
		return 1;

	return 0;
}

int vcn_sparse_solve_CG_precond_Cholesky
(const vcn_sparse_t *const A,
 const double *const b, 
 double *_x,                /* Out */
 uint32_t ktrunc,
 uint32_t max_iter, double tolerance,
 uint32_t* niter_performed,     /* Out (NULL if not required) */
 double* tolerance_reached, /* Out (NULL if not required) */
 uint32_t omp_parallel_threads)
/* Return the num of iterations */
{
	/* Solve Ax = b with Conjugate Gradient preconditioned with 
	 * Cholesky truncated.
	 */
	double* g = (double*)calloc(A->N, sizeof(double));
	double* p = (double*)calloc(A->N, sizeof(double));
	double* q = (double*)calloc(A->N, sizeof(double));
	double* w = (double*)calloc(A->N, sizeof(double));

	double dot_gg = 0;
  
#pragma omp parallel for reduction(+:dot_gg) num_threads(omp_parallel_threads) schedule(guided)
	for(uint32_t i=0; i< A->N; i++){
		double sum = 0;
		for(uint32_t j=0; j< A->rows_size[i]; j++)
			sum += A->rows_values[i][j] * _x[A->rows_index[i][j]];
		g[i] = sum - b[i];
		dot_gg += g[i]*g[i];
	}
	/* Solve H Ht q = g for q */
	vcn_sparse_t *H = sparse_allocate(A->N);
	vcn_sparse_t *Ht = sparse_allocate(A->N);
	cholesky_symbolic(A, H, Ht, ktrunc);

	vcn_sparse_decompose_Cholesky(A,H,Ht, omp_parallel_threads);
	vcn_sparse_solve_LU(H,Ht,g,q);

#pragma omp parallel for num_threads(omp_parallel_threads)
	for(uint32_t i=0; i< A->N; i++)
		p[i] = -q[i];

	uint32_t k = 0;
	while(dot_gg > tolerance*tolerance && k < max_iter){
		double dot_pw = 0;
		double dot_gq = 0;
		dot_gg = 0;
#pragma omp parallel for reduction(+:dot_pw, dot_gg, dot_gq) num_threads(omp_parallel_threads)
		for(uint32_t i=0; i < A->N; i++){
			w[i] = 0;
			for(uint32_t j=0; j< A->rows_size[i]; j++)
				w[i] += A->rows_values[i][j] * p[A->rows_index[i][j]];
			dot_pw += p[i]*w[i];
			dot_gg += g[i]*g[i];
			dot_gq += g[i]*q[i];
		}
		double alphak = dot_gq/dot_pw;
#pragma omp parallel for num_threads(omp_parallel_threads)
		for(uint32_t i=0; i < A->N; i++){
			_x[i] += alphak*p[i];
			g[i] += alphak*w[i];
		}
		/* Solve H Ht q = g for q */
		vcn_sparse_solve_LU(H,Ht,g,q);

		double dot_gkqk = 0;
#pragma omp parallel for reduction(+:dot_gkqk) num_threads(omp_parallel_threads)
		for(uint32_t i=0; i< A->N; i++)
			dot_gkqk += g[i]*q[i];

		double betak = dot_gkqk/dot_gq;
#pragma omp parallel for num_threads(omp_parallel_threads)
		for(uint32_t i=0; i< A->N; i++)
			p[i] = -q[i] + betak*p[i];
		k++;
	}
	/* Free memory */
	free(g);
	free(p);
	free(q);
	free(w);
	vcn_sparse_destroy(H);
	vcn_sparse_destroy(Ht);

	if(niter_performed != NULL) niter_performed[0]= k;

	if(tolerance_reached != NULL) *tolerance_reached = sqrt(dot_gg);

	if(dot_gg > tolerance*tolerance)
		return 1;

	return 0;
}

int vcn_sparse_solve_CG_precond_fsai
(const vcn_sparse_t *const A,
 const double *const b,
 double *_x,                /* Out */
 double threshold,
 uint32_t max_iter,	double tolerance,
 uint32_t* niter_performed,     /* Out (NULL if not required) */
 double* tolerance_reached, /* Out (NULL if not required) */
 uint32_t omp_parallel_threads)
/* Return the num of iterations */
{
	/* Conjugate gradient preconditioned with "Factorized sparse 
	 * approximated inverse" 
	 */
	double *D = (double*)calloc(A->N, sizeof(double));
	double *siD = (double*)calloc(A->N, sizeof(double));

	vcn_sparse_t* G  = sparse_allocate(A->N);
	vcn_sparse_t* Gt = sparse_allocate(A->N);
	/* Generate D diagonal matrix as
	 *           
	 *     Dii = |Aii|,   if |Aii| > 0
	 *            1       otherwise
	 */

#pragma omp parallel for num_threads(omp_parallel_threads)
	for(uint32_t i=0; i < A->N; i++){
		D[i] = vcn_sparse_get(A,i,i);
		if(D[i] == 0)
			D[i] = 1;
		/* Compute D^(-1/2) */
		siD[i] = 1/sqrt(D[i]);
	}
  
	/* Generate structure of G lower triangular matrix 
	 *
	 *    G = 1 ,  if (i == j) or (|[D^(-1/2) A D^(-1/2)]_ij| > threshold)
	 *        0 ,  otherwise
	 *
	 */
	for(uint32_t i=0; i < A->N; i++){
		uint32_t isize = 0;
		uint32_t isizet = 0;

#pragma omp parallel for reduction(+:isize,isizet) num_threads(omp_parallel_threads)
		for(uint32_t q=0; q< A->rows_size[i]; q++){
			uint32_t j = A->rows_index[i][q];
			if(i == j || 
			   fabs(siD[i] * A->rows_values[i][q] * siD[j])
			   > threshold){
				if(i > j)
					isize++;
				else if(i < j)
					isizet++;
				else{
					isize++;
					isizet++;
				}
			}
		}
		G->rows_size[i] = isize;
		G->rows_index[i] = (uint32_t*)calloc(isize,sizeof(uint32_t));
		G->rows_values[i] = (double*)calloc(isize,sizeof(double));

		Gt->rows_size[i] = isizet;
		Gt->rows_index[i] = (uint32_t*)calloc(isizet,sizeof(uint32_t));
		Gt->rows_values[i] = (double*)calloc(isizet,sizeof(double));
	}

#pragma omp parallel for num_threads(omp_parallel_threads)
	for(uint32_t i=0; i < A->N; i++){
		/* Compute values of ~G */
		double* subA = (double*)
			calloc(G->rows_size[i]*G->rows_size[i], sizeof(double));
		/* The data of vector g is not allocated, is a pointer to each row of ~G */
		double* subg = G->rows_values[i];
		double *delta = (double*)calloc(G->rows_size[i], sizeof(double));
		uint32_t k = 0;
		for(uint32_t q = 0; q < A->rows_size[i]; q++){
			uint32_t j = A->rows_index[i][q];
			if(i == j || 
			   fabs(siD[i] * A->rows_values[i][q] * siD[j]) > threshold){
				if(i >= j){
					G->rows_index[i][k] = j;
					for(uint32_t l=0; l<k; l++){
						subA[k*G->rows_size[i] + l] = 
							vcn_sparse_get(A,j,G->rows_index[i][l]);
						subA[l*G->rows_size[i] + k] = 
							vcn_sparse_get(A,G->rows_index[i][l],j);
					}
					subA[k*G->rows_size[i] + k] = vcn_sparse_get(A,j,j);
					if(i == j)
						delta[k] = 1;
					else
						delta[k] = 0;
					k++;
				}
			}
		}
		double* L = (double*) 
			calloc(G->rows_size[i]*G->rows_size[i], sizeof(double));
		vcn_matrix_cholesky_decomposition(subA, L, G->rows_size[i]);
		vcn_matrix_cholesky_solve(L, delta, subg, G->rows_size[i]);
		/* Finally do G = [~G]*D   */
		for(uint32_t q=0; q < G->rows_size[i]; q++)
			G->rows_values[i][q] *= D[G->rows_index[i][q]];

		/* Free memory */
		free(subA);
		free(L);
		free(delta);
	}
	/* Store G transposed */
	vcn_sparse_transpose(G,Gt);

	/* Free memory */
	free(D);
	free(siD);

	/* Solve Ax = b with Conjugate Gradient method */
	double* r = (double*)calloc(A->N, sizeof(double));
	double* p = (double*)calloc(A->N, sizeof(double));
	double* w = (double*)calloc(A->N, sizeof(double));
	double* Gr = (double*)calloc(A->N, sizeof(double));
	double* Mr = (double*)calloc(A->N, sizeof(double));

	double dot_rr = 0;

#pragma omp parallel for reduction(+:dot_rr) schedule(guided) num_threads(omp_parallel_threads)
	for(uint32_t i=0; i< A->N; i++){
		double Ax_i = 0;
		for(uint32_t j=0; j< A->rows_size[i]; j++)
			Ax_i += A->rows_values[i][j] * _x[A->rows_index[i][j]];
		r[i] = b[i] - Ax_i;
		dot_rr += r[i]*r[i];
	}
	/* Compute Gr */

#pragma omp parallel for num_threads(omp_parallel_threads)
	for(uint32_t i=0; i< A->N; i++){
		Gr[i] = 0;
		for(uint32_t j=0; j< G->rows_size[i]; j++)
			Gr[i] += G->rows_values[i][j] * r[G->rows_index[i][j]];
	}
	/* Compute Mr <- G'(Gr) */

#pragma omp parallel for num_threads(omp_parallel_threads)
	for(uint32_t i=0; i< A->N; i++){
		Mr[i] = 0;
		for(uint32_t j=0; j< Gt->rows_size[i]; j++)
			Mr[i] += Gt->rows_values[i][j] * Gr[Gt->rows_index[i][j]];
		p[i] = Mr[i];
	}
	uint32_t k = 0;
	/* Start iterations */
	while(dot_rr > tolerance*tolerance && k < max_iter){
		double dot_pw = 0;
		double dot_rMr = 0;

#pragma omp parallel for reduction(+:dot_pw, dot_rMr) num_threads(omp_parallel_threads) schedule(guided)
		for(uint32_t i=0; i< A->N; i++){
			w[i] = 0;
			for(uint32_t j=0; j< A->rows_size[i]; j++)
				w[i] += A->rows_values[i][j] * p[A->rows_index[i][j]];
			dot_pw += p[i]*w[i];
			dot_rMr += r[i]*Mr[i];
		}
		double alphak = dot_rMr/dot_pw;
		dot_rr = 0;

#pragma omp parallel for reduction(+:dot_rr) num_threads(omp_parallel_threads) schedule(guided)
		for(uint32_t i=0; i< A->N; i++){
			_x[i] += alphak*p[i];
			r[i] -= alphak*w[i];
			dot_rr += r[i]*r[i];
		}
		/* Compute Gr */

#pragma omp parallel for num_threads(omp_parallel_threads)
		for(uint32_t i=0; i< A->N; i++){
			Gr[i] = 0;
			for(uint32_t j=0; j< G->rows_size[i]; j++)
				Gr[i] += G->rows_values[i][j] * r[G->rows_index[i][j]];
		}
		/* Compute Mr <- G'(Gr) */
		double dot_rkMrk = 0;

#pragma omp parallel for reduction(+:dot_rkMrk) num_threads(omp_parallel_threads)
		for(uint32_t i=0; i< A->N; i++){
			Mr[i] = 0;
			for(uint32_t j=0; j< Gt->rows_size[i]; j++)
				Mr[i] += Gt->rows_values[i][j] * Gr[Gt->rows_index[i][j]];
			dot_rkMrk += r[i]*Mr[i];
		}
		double betak = dot_rkMrk/dot_rMr;

#pragma omp parallel for num_threads(omp_parallel_threads)
		for(uint32_t i=0; i< A->N; i++)
			p[i] = Mr[i] + betak*p[i];
		k++;
	}
	/* Free memory */
	vcn_sparse_destroy(G);
	vcn_sparse_destroy(Gt);
	free(r);
	free(p);
	free(w);
	free(Gr);
	free(Mr);

	if(niter_performed != NULL) niter_performed[0]= k;

	if(tolerance_reached != NULL) *tolerance_reached = sqrt(dot_rr);

	if(dot_rr > tolerance*tolerance)
		return 1;

	return 0;
}

typedef struct{
	/* OPPORTUNITY: Use libre_dstructs */
	/* Symbolic Cholesky Set
	 * This structure is going to be useful to create
	 * lists that represent sets to implement the 
	 * "symbolic cholesky factorization" algorithm.
	 */
	uint32_t index;
	void *next;
}sc_set;

static void cholesky_symbolic(const vcn_sparse_t * const A, 
			      vcn_sparse_t *_L, vcn_sparse_t* _Lt, 
			      uint32_t ktrunc)
{
	/* The algorithm is going to store in "l" the index of
	 * columns that must have each column of the lower triangular
	 * matrix "L" from A = LL' (without the index of the column "Lii"), 
	 * using only 'ktrunc' diagonals from A.
	 * It could look something like:
	 *   l[0]:   1 -> 2 -> 3 -> 5 -> NULL
	 *   l[1]:   2 -> 4 -> 6 -> 7 -> 8 -> NULL
	 *   l[2]:   3 -> 5 -> 7 -> 9 -> NULL
	 *     :
	 *     :
	 *   l[n-2]: n-1 -> NULL
	 * for a matrix of size n, starting the index from 0.
	 * After each iteration, the algortihm is going to allocate the
	 * i row of the upper matrix "Lt", finally, the matrix "L" will
	 * be allocated in row compress format (vcn_sparse_t structure).
	 *
	 * Reference:
	 *  GALLIVAN, Heath, Ng, Ortega, Peyton, Plemmons, Romine, 
	 *    Sameh & Voigt. 
	 *    "Parallel Algorithms for Matrix Computations"
	 *    SIAM 1990 p86-88
	 */
	uint32_t *L_size = calloc(A->N, sizeof(*L_size));

	sc_set** r = calloc(A->N, sizeof(*r));

	for (uint32_t j = 0; j < A->N; j++) { 
		sc_set* lj = NULL;
		uint32_t lj_size = 1;
		uint32_t _i = sparse_bsearch_row(A,j,j, 0, A->rows_size[j]-1);

		/* lj <- aj ************************************************/
		sc_set *iterator_lj;                                     /**/
		for(uint32_t i = _i+1; i<A->rows_size[j]; i++){     /**/
			sc_set* aji = (sc_set*)malloc(sizeof(sc_set));         /**/
			aji->index = A->rows_index[j][i];                      /**/
			aji->next = NULL;                                      /**/
			if(lj == NULL){                                        /**/
				lj = aji;                                            /**/
				iterator_lj = lj;                                    /**/
			}else{                                                 /**/
				iterator_lj->next = aji;                             /**/
				iterator_lj = aji;                                   /**/
			}                                                      /**/
			lj_size++;                                             /**/
		}                                                        /**/
		/***********************************************************/

		sc_set* iterator_rj = r[j];
		while(iterator_rj != NULL){
			/* lj <- lj U {x in li|x-j<k} \ {j} ***********************/
			/*                         ^ This k refers to ktrunc     **/
			/*                           from parameters.            **/
			uint32_t i = iterator_rj->index;                            /**/
			if(lj != NULL){                                         /**/
				for(uint32_t k=1; k<_Lt->rows_size[i]; k++){     /**/
					uint32_t index = _Lt->rows_index[i][k];                 /**/
					if(index != j && index - j < ktrunc){               /**/
						sc_set* next_itr_lj = lj;                         /**/
						sc_set* iterator_lj = NULL;                       /**/
						short flag1 = 0;          /* {index} > all in lj   **/
						while(next_itr_lj != NULL && flag1 == 0){         /**/
							if(next_itr_lj->index > index)                  /**/
								flag1 = 1;            /* {index} < some in lj  **/
							else if(next_itr_lj->index == index)            /**/
								flag1 = 2;            /* {index} in j          **/
							else{                                           /**/
								iterator_lj = next_itr_lj;                    /**/
								next_itr_lj = (sc_set*)next_itr_lj->next;     /**/
							}                                               /**/
						}                                                 /**/
						if(flag1 != 2){                                   /**/
							sc_set* node = (sc_set*)malloc(sizeof(sc_set)); /**/
							node->index = index;                            /**/
							node->next = next_itr_lj;                       /**/
							if(iterator_lj != NULL)                         /**/
								iterator_lj->next = node;                     /**/
							else                                            /**/
								lj = node;                                    /**/
							lj_size++;                                      /**/
						} /* end if (flag1 != 2)                           **/
					}                                                   /**/
				}                                                     /**/
			}else{                                                  /**/
				short flag = 0;             /* {index} == {j}          **/
				uint32_t k = 1;                                           /**/
				while(k < _Lt->rows_size[i] && flag == 0){            /**/
					if(_Lt->rows_index[i][k] != j)                      /**/
						flag = 1;               /* {index} != {j}          **/
					else                                                /**/
						k++;                                              /**/
				}                                                     /**/
				if(flag == 1){                                        /**/
					lj = (sc_set*)malloc(sizeof(sc_set));               /**/
					lj->index = _Lt->rows_index[i][k];                  /**/
					lj->next = NULL;                                    /**/
					lj_size++;                                          /**/
					k++;                                                /**/
				}                                                     /**/
				sc_set* iterator_lj = lj;                             /**/
				while(k < _Lt->rows_size[i]){                         /**/
					if(_Lt->rows_index[i][k] != j){                     /**/
						sc_set* node = (sc_set*)malloc(sizeof(sc_set));   /**/
						node->index = _Lt->rows_index[i][k];              /**/
						node->next = NULL;                                /**/
						iterator_lj->next = node;                         /**/
						lj_size++;                                        /**/
						iterator_lj = (sc_set*)iterator_lj->next;         /**/
					}                                                   /**/
					k++;                                                /**/
				}                                                     /**/
			}                                                       /**/
			/**********************************************************/
			iterator_rj = (sc_set*)iterator_rj->next;
		}

		if(lj != NULL){
			/* p <- min(i in lj) if (lj != NULL) */
			uint32_t p = lj->index;

			/* rp <- rp U {j} *****************************************/
			if(r[p] != NULL){                                       /**/
				sc_set* iterator_rp = r[p];                           /**/
				short flag = 0;                /* {j} is not in rp     **/
				while(iterator_rp->next != NULL && flag == 0){        /**/
					if(iterator_rp->index == j)                         /**/
						flag = 1;                  /* {i} in rp            **/
					iterator_rp = (sc_set*)iterator_rp->next;           /**/
				}                                                     /**/
				if(flag == 0){                                        /**/
					if(iterator_rp->index != j){                        /**/
						sc_set* node = (sc_set*)malloc(sizeof(sc_set));   /**/
						node->index = j;                                  /**/
						node->next = NULL;                                /**/
						iterator_rp->next = node;                         /**/
					}                                                   /**/
				}                                                     /**/
			}else{                                                  /**/
				r[p] = (sc_set*)malloc(sizeof(sc_set));               /**/
				r[p]->index = j;                                      /**/
				r[p]->next = NULL;                                    /**/
			}                                                       /**/
			/**********************************************************/
		}

		/*************** Allocate the jth row of "Lt" **********************/
		_Lt->rows_values[j] = (double*)calloc(lj_size,sizeof(double));   /**/
		_Lt->rows_index[j] = (uint32_t*)calloc(lj_size,sizeof(uint32_t));        /**/
		_Lt->rows_size[j] = lj_size;                                     /**/
		iterator_lj = lj;                                                /**/
		_Lt->rows_index[j][0] = j;                                       /**/
		for(uint32_t k=1; k<lj_size; k++){                          /**/
			_Lt->rows_index[j][k] = iterator_lj->index;                    /**/
			L_size[iterator_lj->index] += 1;                               /**/
			iterator_lj = (sc_set*)iterator_lj->next;                      /**/
		}                                                                /**/
		/* Free memory */                                                /**/
		while(lj != NULL){                                               /**/
			sc_set* lj_free = lj;                                          /**/
			lj = (sc_set*)lj->next;                                        /**/
			free(lj_free);                                                 /**/
		}                                                                /**/
		/*******************************************************************/
    
	}
	/****************************** Allocate "L" ****************************/
	for(uint32_t i=0; i<_L->N; i++){                                 /**/
		L_size[i]++;            /* To include main diagonal */              /**/
		_L->rows_values[i] = (double*)calloc(L_size[i],sizeof(double));     /**/
		_L->rows_index[i] = (uint32_t*)calloc(L_size[i],sizeof(uint32_t));          /**/
		_L->rows_size[i] = L_size[i];                                       /**/
	}                                                                     /**/
	/* Cycle to set columns in L */                                       /**/
	uint32_t* L_index = (uint32_t*)calloc(_L->N,sizeof(uint32_t));                    /**/
	/* This "for" can not be parallelizad to guaranty increasing sort */  /**/
	for(uint32_t i=0; i<_Lt->N; i++){                                /**/
		for(uint32_t j=0; j<_Lt->rows_size[i]; j++){                   /**/
			uint32_t index = _Lt->rows_index[i][j];                               /**/
			_L->rows_index[index][L_index[index]++] = i;                      /**/
		}                                                                   /**/
	}                                                                     /**/
	/************************************************************************/

	/* Free memory */
	free(L_size);
	free(L_index);
	for(uint32_t i=0; i< A->N; i++){
		sc_set* iterator_ri = r[i];
		while(iterator_ri != NULL){
			sc_set* rm = iterator_ri;
			iterator_ri = (sc_set*)iterator_ri->next;
			free(rm);
		}
	}
	free(r);
}


int vcn_sparse_alloc_LU
(const vcn_sparse_t *const restrict A,
 vcn_sparse_t** L, vcn_sparse_t** U)
{
	*L = sparse_allocate(A->N);
	*U = sparse_allocate(A->N);
  
	/* Run cholesky symbolic to allocate L and Lt */
	cholesky_symbolic(A, *L, *U, A->N);

	return 0;
}

int vcn_sparse_decompose_Cholesky(const vcn_sparse_t *const A,
				  vcn_sparse_t *L,             /* Out */
				  vcn_sparse_t* Lt,            /* Out */
				  uint32_t omp_parallel_threads)
{
	/* "L" must be a lower triangular matrix with the main diagonal 
	 * complete, and "Lt" must be an upper triangular matrix with
	 * the main diagonal complete.
	 * The structure of L must be congrous with Lt, since Lt = L'.
	 */
    
	/* Compute the decomposition */
	for(uint32_t j=0; j< A->N; j++){
		L->rows_values[j][L->rows_size[j]-1] = vcn_sparse_get(A, j, j);
    
		double sum = 0;
		for(uint32_t q = 0; q < L->rows_size[j]-1; q++)
			sum += POW2(L->rows_values[j][q]);
      
		L->rows_values[j][L->rows_size[j]-1] -= sum;

		if(L->rows_values[j][L->rows_size[j]-1] <= 0.0)
			return 1;
              
		double valuejj = sqrt(L->rows_values[j][L->rows_size[j]-1]);
		L->rows_values[j][L->rows_size[j]-1] = valuejj;
		Lt->rows_values[j][0] = valuejj;

#pragma omp parallel for  num_threads(omp_parallel_threads) schedule(guided)
		for(uint32_t q = 1; q < Lt->rows_size[j]; q++){
			uint32_t i = Lt->rows_index[j][q];
			/*** L_ij <- A_ij ********************************************************/
			uint32_t L_jindex = sparse_bsearch_row(L, i, j, 0, L->rows_size[i]-1);     /**/
			L->rows_values[i][L_jindex] = vcn_sparse_get(A, i, j);                     /**/
			/*************************************************************************/
			uint32_t r = 0;
			uint32_t s = 0;
			uint32_t _ro = L->rows_index[i][r];
			uint32_t _sigma = L->rows_index[j][s];
			bool flag = true;   /* Flag to know when to stop the cylce */
			while(flag){
				while(_ro < _sigma)
					_ro = L->rows_index[i][++r];
				while(_ro > _sigma)
					_sigma = L->rows_index[j][++s];
				while(_ro == _sigma){
					if(_ro == j){
						flag = false;   /* Finish the cycle */
						break;
					}
					double vir = L->rows_values[i][r];
					double vjs = L->rows_values[j][s];
					L->rows_values[i][L_jindex] -= vir*vjs;
					_ro = L->rows_index[i][++r];
					_sigma = L->rows_index[j][++s];
				}
			}

			L->rows_values[i][L_jindex] /= L->rows_values[j][L->rows_size[j]-1];
			Lt->rows_values[j][q] = L->rows_values[i][L_jindex];
		}
	}

	/* Successful exit */
	return 0;
}

void vcn_sparse_decompose_LU(const vcn_sparse_t *const Ar,
			     vcn_sparse_t *L, vcn_sparse_t* U,
			     uint32_t omp_parallel_threads){

	/* Create Ut to compute faster the decomposition */
	vcn_sparse_t* Ut = vcn_sparse_clone(L);

	/* Compute the decomposition */
	for(uint32_t j=0; j< Ar->N; j++){
		L->rows_values[j][L->rows_size[j]-1] = 1.0;
		U->rows_values[j][0] = vcn_sparse_get(Ar, j, j);

		double sum = 0;

#pragma omp parallel for schedule(guided) reduction(+:sum) num_threads(omp_parallel_threads)
		for(uint32_t q=0; q< L->rows_size[j]-1; q++)
			sum += L->rows_values[j][q] * Ut->rows_values[j][q];
    

		U->rows_values[j][0] -= sum;
		Ut->rows_values[j][Ut->rows_size[j]-1] = U->rows_values[j][0];
        
#pragma omp parallel for schedule(guided) num_threads(omp_parallel_threads)
		for(uint32_t q = 1; q < U->rows_size[j]; q++){
			uint32_t i = U->rows_index[j][q];
			/*** L_ij <- A_ij *******************************************************/
			uint32_t L_jindex = sparse_bsearch_row(L, i, j, 0, L->rows_size[i]-1);/**/
			L->rows_values[i][L_jindex] = vcn_sparse_get(Ar, i, j);               /**/
			/************************************************************************/
			/*** U_ji <- A_ji *******************************************************/
			U->rows_values[j][q] = vcn_sparse_get(Ar, j, i);                      /**/
			/************************************************************************/
			register uint32_t r = 0;
			register uint32_t s = 0;
			register uint32_t _ro =    L->rows_index[i][r];
			register uint32_t _sigma = L->rows_index[j][s];
			bool flag = true;   /* Flag to know when to stop the cylce */
			while(flag){
				while(_ro < _sigma)
					_ro = L->rows_index[i][++r];
				while(_ro > _sigma)
					_sigma = L->rows_index[j][++s];
				while(_ro == _sigma){
					if(_ro == j){
						flag = false;   /* Finish the cycle */
						break;
					}
					double vir = L->rows_values[i][r];
					double vjs = Ut->rows_values[j][s];
					L->rows_values[i][L_jindex] -= vir*vjs;

					vjs = L->rows_values[j][s];
					vir = Ut->rows_values[i][r];
					U->rows_values[j][q] -= vir*vjs;

					_ro = L->rows_index[i][++r];
					_sigma = L->rows_index[j][++s];
				}
			}

			L->rows_values[i][L_jindex] /= U->rows_values[j][0];
      
			Ut->rows_values[i][L_jindex] = U->rows_values[j][q];
		}
	}
	/* Free memory */
	vcn_sparse_destroy(Ut);
}

void  vcn_sparse_solve_LU(const vcn_sparse_t *const L, 
			  const vcn_sparse_t *const U,
			  const double *const b,
			  double* _x  /* Out */)
{
	double* z = calloc(L->N, sizeof(*z));
	vcn_sparse_forward_solve(L, b, z);
	vcn_sparse_backward_solve(U, z, _x);
	free(z);
}

int vcn_sparse_solve_Cholesky(const vcn_sparse_t *const A,
			      const double *const b,
			      double* x,  /* Out */
			      uint32_t omp_parallel_threads){
	vcn_sparse_t *L = NULL; 
	vcn_sparse_t *U = NULL;
	vcn_sparse_alloc_LU(A, &L, &U);
	if (NULL == L)
		return 10;
	
	int solver_status = vcn_sparse_decompose_Cholesky(A, L, U, 
							  omp_parallel_threads);

	if (0 == solver_status)
		vcn_sparse_solve_LU(L, U, b, x);

	vcn_sparse_destroy(L);
	vcn_sparse_destroy(U);

	return solver_status;
}

int vcn_sparse_solve_using_LU(const vcn_sparse_t *const A,
			      const double *const b,
			      double* x,  /* Out */
			      uint32_t omp_parallel_threads)
{
	vcn_sparse_t *L = NULL; 
	vcn_sparse_t *U = NULL;
	vcn_sparse_alloc_LU(A, &L, &U);
	if(NULL == L)
		return 1;

	vcn_sparse_decompose_LU(A, L, U, omp_parallel_threads);

	vcn_sparse_solve_LU(L, U, b, x);

	vcn_sparse_destroy(L);
	vcn_sparse_destroy(U);

	return 0;
}

void vcn_sparse_forward_solve(const vcn_sparse_t *const L,
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
			_x[i] /= vcn_sparse_get(L, i, i);
	}
}

void vcn_sparse_backward_solve(const vcn_sparse_t *const U,
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
			idx = sparse_bsearch_row(U, i, i, 0, U->rows_size[i]-1);
		for(int q = idx+1; q < U->rows_size[i]; q++)
			_x[i] -= U->rows_values[i][q]*_x[U->rows_index[i][q]];
    
		_x[i] /= U->rows_values[i][idx];
	}
}

void vcn_matrix_2X2_inverse(const double *const A,
			    double* A_inv)
{
	double det = (A[0] * A[3] - A[1] * A[2]);
	A_inv[0] = A[3] / det;
	A_inv[1] = -A[1] / det;
	A_inv[2] = -A[2] / det;
	A_inv[3] = A[0] / det;
}

void vcn_matrix_2X2_eigen(const double *const A,
			  double* Lambda, /* Output (diag) */
			  double* P,      /* Output */
			  double tolerance)
/* Returns A decomposed into P Lambda P' */
{
	if (fabs(A[1]) < tolerance && fabs(A[2]) < tolerance) {
		Lambda[0] = A[0];
		Lambda[1] = A[3];
		P[0] = 1.0;
		P[2] = 0.0;
		P[1] = 0.0;
		P[3] = 1.0;
	} else {
		double T = A[0] + A[3];
		double D = A[0] * A[3] - A[1] * A[2];
		double root = sqrt(POW2(T)/4.0 - D);
		Lambda[0] = T/2.0 + root;
		Lambda[1] = T/2.0 - root;
		if (fabs(A[2]) > tolerance) {
			P[0] = Lambda[0] - A[3];
			P[2] = A[2];
			P[1] = Lambda[1] - A[3];
			P[3] = A[2];
		} else {
			P[0] = A[1];
			P[2] = Lambda[0] - A[0];
			P[1] = A[1];
			P[3] = Lambda[1] - A[0];
		}
	}
	if (fabs(Lambda[0]) < fabs(Lambda[1])) {
		double aux = Lambda[0];
		Lambda[0] = Lambda[1];
		Lambda[1] = aux;

		aux = P[0];
		P[0] = P[1];
		P[1] = aux;

		aux = P[2];
		P[2] = P[3];
		P[3] = aux;
	}
	/* Normalize eigenvectors */
	double normalizer = sqrt(POW2(P[0]) + POW2(P[2]));
	P[0] /= normalizer;
	P[2] /= normalizer;
	normalizer = sqrt(POW2(P[1]) + POW2(P[3]));
	P[1] /= normalizer;
	P[3] /= normalizer;
}
double vcn_matrix_2X2_det(double *A){
	return A[0] * A[3] - A[1] * A[2];
}

double vcn_matrix_3X3_det(double *A){
	return 
		A[0] * A[4] * A[8] + 
		A[3] * A[7] * A[2] +
		A[6] * A[1] * A[5] -
		A[2] * A[4] * A[6] -
		A[5] * A[7] * A[0] -
		A[8] * A[1] * A[3];
}

void vcn_matrix_2X2_inverse_destructive(double *A){
	double det = A[0] * A[3] - A[1] * A[2];
	double tmp = A[0];
	A[0] = A[3]/det;
	A[3] = tmp/det;
	A[1] = -A[1]/det;
	A[2] = -A[2]/det;
}

void vcn_matrix_3X3_inverse_destructive(double *A){
	double det = vcn_matrix_3X3_det(A);
	double a11, a12, a13, a21, a22, a23, a31, a32, a33;
	a11 = A[0];
	a12 = A[1];
	a13 = A[2];
	a21 = A[3];
	a22 = A[4];
	a23 = A[5];
	a31 = A[6];
	a32 = A[7];
	a33 = A[8];
	A[0] = (a33*a22-a32*a23)/det;
	A[1] = -(a33*a12-a32*a13)/det;
	A[2] = (a23*a12-a22*a13)/det;
	A[3] = -(a33*a21-a31*a23)/det;
	A[4] = (a33*a11-a31*a13)/det;
	A[5] = -(a23*a11-a21*a13)/det;
	A[6] = (a32*a21-a31*a22)/det;
	A[7] = -(a32*a11*a31*a12)/det;
	A[8] = (a22*a11-a21*a12)/det;
}

int vcn_matrix_cholesky_decomposition
(const double *const A,
 double* _LplusLt,      /* Out */
 uint32_t N){
	/* The 'L' computed correspond to [L + Lt] */
	for (uint32_t j=0; j < N; j++) {
		_LplusLt[j*N+j] = A[j*N+j];
		for (uint32_t k=0; k<j; k++)
			_LplusLt[j*N+j] -= _LplusLt[j*N+k]*_LplusLt[j*N+k];
		if (_LplusLt[j*N+j] <= 1e-16)
			return 1;
		_LplusLt[j*N+j] = sqrt(_LplusLt[j*N+j]);
		for (uint32_t i=j+1;  i<N; i++) {
			_LplusLt[i*N+j] = A[i*N+j];
			for (uint32_t k=0; k<j; k++)
				_LplusLt[j*N+j] -= _LplusLt[i*N+k]*_LplusLt[j*N+k];
			_LplusLt[i*N+j] /= _LplusLt[j*N+j];
			_LplusLt[j*N+i] = _LplusLt[i*N+j];
		}
	}
	return 0;
}

void vcn_matrix_cholesky_solve
(const double *const LplusLt,
 const double *const b, 
 double* _x,            /* Out */
 uint32_t N)
/* Solve the system LL'x = b, where LL'= A */
{
	double* z = (double*)calloc(N, sizeof(double));
	vcn_matrix_forward_solve(LplusLt, b, z, N);
	vcn_matrix_backward_solve(LplusLt, z, _x, N);
	/* Free memory */
	free(z);
}

double vcn_matrix_cond1(const double *const A, int n){
	double *Acopy = (double*)malloc(n*n*sizeof(double));
	memcpy(Acopy, A, n*n*sizeof(double));
	double *x = (double*)malloc(n*sizeof(double));
	double *p = (double*)malloc(n*sizeof(double));
	double *pm = (double*)malloc(n*sizeof(double));
	/* Compute QR decomposition */
	int sing;
	double *c = (double*)malloc(n*sizeof(double));
	double *diag = (double*)malloc(n*sizeof(double));
	vcn_matrix_qr_decomposition(Acopy, n, c, diag, &sing);
	/* Compute ||A||_1 (Max absolute column sum) */
	double estimation = fabs(diag[0]);
	for(uint32_t i=1; i<n; i++){
		double sum = 0;
		for(uint32_t j=0; j<i; j++)
			sum +=  fabs(Acopy[i*n+j]);
		double tmp = fabs(diag[i]) + sum;
		estimation = (tmp>estimation)?tmp:estimation;
	}
	/* Solve R^T x = e, selecting e as they proceed */
	x[0] = 1.0/diag[0];
	for(uint32_t i=1; i<n; i++)
		p[i] = Acopy[i]*x[0];
	for(uint32_t i=1; i<n; i++){
		/* Select ej and calculate xj */
		double xp = (1-p[i])/diag[i];
		double xm = (-1-p[i])/diag[i];
		double tmp = fabs(xp);
		double tmpm = fabs(xm);
		for(uint32_t j=i+1; j<n; j++){
			pm[j] = p[j] + Acopy[i*n+j]*xm;
			tmpm += (fabs(pm[j])/fabs(diag[j]));
			p[j] = p[j] + Acopy[i*n+j]*xp;
			tmp += (fabs(p[j])/fabs(diag[j]));
		}
		if(tmp > tmpm) x[i] = xp;    /* ej = 1 */
		else{
			/* ej = -1 */
			x[i] = xm;
			for(uint32_t j=i+1; j<n; j++)
				p[j]  = pm[j];
		}
	}
	double xnorm = 0;
	for(uint32_t i=0; i<n; i++) xnorm += fabs(x[i]);
	estimation /= xnorm;
	/* Solve Ry = x */
	vcn_matrix_qr_solve(Acopy, n, c, diag, x);

	xnorm = 0;
	for(uint32_t i=0; i<n; i++) xnorm += fabs(x[i]);
	estimation *= xnorm;
	/* Free memory */
	free(Acopy);
	free(x);
	free(p);
	free(pm);
	free(diag);
	free(c);
	return estimation;
}

double vcn_matrix_cond2(const double *const A, int n){
	double *Acopy = (double*)malloc(n*n*sizeof(double));
	memcpy(Acopy, A, n*n*sizeof(double));
	double* w = (double*)malloc(n*sizeof(double));
	double* V = (double*)malloc(n*n*sizeof(double));
	/* Compute SVD Decomposition */ 
	vcn_matrix_svd_decomposition(Acopy, w, V, n, n);
	/* Compute condition number estimation */
	double estimation = w[0]/w[n-1];
	/* Free memory */
	free(Acopy);
	free(w);
	free(V);
	return estimation;
}

void vcn_matrix_qr_decomposition(double *A, /* Overwritten */
				 int n,
				 double *c, double *d, int *sing){
	/* Numerical Recipes in C.
	 * Constructs the QR decomposition of A. The upper triangular matrix
	 * R is returned in the upper triangle of A, except for the diagonal 
	 * elements of R which are returned in d[1..n]. The orthogonal matrix 
	 * Q is represented as a product of n-1 Householder matrices Q_1...Q_{n-1},
	 * where Q_j = 1 - u_j prod u_j/c_j. The ith component of u_j is zero
	 * for i=1,...,j-1 while the nonzero components are returned in A_ij 
	 * for i=j,...,n.  "sing" return as true(1) if singularity is encountered
	 * during the decomposition, but the decomposition is still completed in 
	 * this case, otherwise it returns false(0).
	 */
	double scale, sigma, sum, tau;
	*sing = 0;
	for(int k=0; k<n-1; k++){
		scale = 0.0;
		for(int i=k; i<n; i++) 
			scale = (scale > fabs(A[i*n+k]))?(scale):(fabs(A[i*n+k]));
		if(scale == 0.0){
			*sing = 1;
			c[k] = d[k] = 0.0;
		}else{
			for(int i=k; i<n; i++)
				A[i*n+k] /= scale;
			sum = 0.0;
			for(int i=k; i<n; i++) 
				sum += A[i*n+k]*A[i*n+k];
			sigma = (A[k*n+k] >= 0.0 ? fabs(sqrt(sum)) : -fabs(sqrt(sum)));
			A[k*n+k] += sigma;
			c[k] = sigma*A[k*n+k];
			d[k] = -scale*sigma;
			for(int j=k+1; j<n; j++){
				sum = 0.0;
				for(int i=k; i<n; i++)
					sum += A[i*n+k]*A[i*n+j];
				tau = sum/c[k];
				for(int i=k; i<n; i++)
					A[i*n+j] -= tau*A[i*n+k];
			}
		}
	}
	d[n-1] = A[n*n-1];
	if(d[n-1] == 0.0) *sing = 1;
}

static inline void matrix_r_solver(const double *const A, int n, double *d, double *b){
	b[n-1] /= d[n-1];
	for(int i=n-2; i>=0; i--){
		double sum = 0.0;
		for(int j=i+1; j<n; j++) sum += A[i*n+j]*b[j];
		b[i] = (b[i]-sum)/d[i];
	}
}

void vcn_matrix_qr_solve(const double *const A,
			 int n, double *c, double *d,
			 double *b /* Solution overwritten */){
	/* Numerical Recipes in C
	 * Solves the set of n linear equations Ax = b. A, c and d are the 
	 * input as the output of the routine vcn_matrix_qr_decomposition  and
	 * are not modified. b is input as the right hand side vector, and 
	 * is overwritten with the solution vector on output.
	 */
	for(int j=0; j<n-1; j++){
		double sum = 0.0;
		for(int i=j; i<n; i++)
			sum += A[i*n+j]*b[i];
		double tau = sum/c[j];
		for(int i=j; i<n; i++)
			b[i] -= tau*A[i*n+j];
	}
	matrix_r_solver(A, n, d, b);
}

void vcn_matrix_svd_decomposition(double *A, /* Overwritten wit U */
				  double *w, double *V, 
				  int n, int m){
	/* Numerical Recipes in C.
	 * Given  a  matrix  A,  this  routine  computes its singular value 
	 * decomposition, A = UWV^T. The matrix U replaces A on output. The
	 * diagonal  matrix  of  singular  values W is output as a vector w. 
	 * The matrix V (not the transpose V^T) is output as V.
	 */
	int flag, i, its, j, jj, k, l, nm;
	double anorm, c, f, g, h, s, scale, x, y, z, *rv1;
	rv1 = (double*)calloc(n, sizeof(double));
	g = scale = anorm = 0.0;
	for(i=0; i<n; i++){
		l = i+1;
		rv1[i] = scale*g;
		g = s = scale = 0.0;
		if(i < m){
			for(k=i; k<m; k++) scale += fabs(A[k*n+i]);
			if(scale){
				for(k=i; k<m; k++){
					A[k*n+i] /= scale;
					s += A[k*n+i]*A[k*n+i];
				}
				f = A[i*n+i];
				g = -(f >= 0.0 ? fabs(sqrt(s)) : -fabs(sqrt(s)));
				h = f*g-s;
				A[i*n+i] = f-g;
				for(j=l; j<n; j++){
					for(s=0.0, k=i; k<m; k++) s += A[k*n+i]*A[k*n+j];
					f = s/h;
					for(k=i; k<m; k++) A[k*n+j] += f*A[k*n+i];
				}
				for(k=i; k<m; k++) A[k*n+i] *= scale;
			}
		}
		w[i] = scale*g;
		g = s = scale = 0.0;
		if(i < m && i != n-1){
			for(k=l; k<n; k++) scale += fabs(A[i*n+k]);
			if(scale){
				for(k=l; k<n; k++){
					A[i*n+k] /= scale;
					s += A[i*n+k]*A[i*n+k];
				}
				f = A[i*n+l];
				g = -(f >= 0.0 ? fabs(sqrt(s)) : -fabs(sqrt(s)));
				h = f*g-s;
				A[i*n+l] = f-g;
				for(k=l; k<n; k++) rv1[k] = A[i*n+k]/h;
				for(j=l; j<m; j++){
					for(s = 0.0, k=l; k<n; k++) s += A[j*n+k]*A[i*n+k];
					for(k=l; k<n; k++) A[j*n+k] += s*rv1[k];
				}
				for(k=l; k<n; k++) A[i*n+k] *= scale;
			}
		}
		anorm = 
			(anorm >= (fabs(w[i])+fabs(rv1[i]))) ? anorm: (fabs(w[i])+fabs(rv1[i]));
	}
	for(i=n-1; i>=0; i--){
		if(i<n-1){
			if(g){
				for(j=l; j<n; j++)
					V[j*n+i] = (A[i*n+j]/A[i*n+l])/g;
				for(j=l; j<n; j++){
					for(s=0.0, k=l; k<n; k++) s += A[i*n+k]*V[k*n+j];
					for(k=l; k<n; k++) V[k*n+j] += s*V[k*n+i];
				}
			}
			for(j=l; j<n; j++) V[i*n+j] = V[j*n+i] = 0.0;
		}
		V[i*n+i] = 1.0;
		g = rv1[i];
		l = i;
	}
	for(i=((n<m)?(n-1):(m-1)); i>=0; i--){
		l = i+1;
		g = w[i];
		for(j=l; j<n; j++) A[i*n+j] = 0.0;
		if(g){
			g = 1.0/g;
			for(j=l; j<n; j++){
				for(s=0.0, k=l; k<m; k++) s += A[k*n+i]*A[k*n+j];
				f = (s/A[i*n+i])*g;
				for(k=i; k<m; k++) A[k*n+j] += f*A[k*n+i];
			}
			for(j=i; j<m; j++) A[j*n+i] *= g;
      
		}else for(j=i; j<n; j++) A[j*n+i] = 0.0;
		++A[i*n+i];
	}
	for(k=n-1; k>=0; k--){
		for(its=0; its<30; its++){
			flag = 1;
			for(l=k; l>=0; l--){
				nm = l-1;
				if(fabs(rv1[l])+anorm == anorm){
					flag = 0;
					break;
				}
				if(fabs(w[nm])+anorm == anorm) break;
			}
			if(flag){
				c = 0.0;
				s = 1.0;
				for(i=l; i<k; i++){
					f = s*rv1[i];
					rv1[i] = c*rv1[i];
					if(fabs(f)+anorm == anorm) break;
					g = w[i];
					h = vcn_math_hypo(f,g);
					w[i] = h;
					h = 1.0/h;
					c = g*h;
					s = -f*h;
					for(j=0; j<m; j++){
						y = A[j*n+nm];
						z = A[j*n+i];
						A[j*n+nm] = y*c+z*s;
						A[j*n+i] = z*c-y*s;
					}
				}
			}
			z = w[k];
			if(l == k){
				if(z < 0.0){
					w[k] = -z;
					for(j=0; j<n; j++) 
						V[j*n+k] = -V[j*n+k];
				}
				break;
			}
			if(its == 29){
				printf("SVD Error 1\n");
				exit(1);
			}
			x = w[l];
			nm = k-1;
			y = w[nm];
			g = rv1[nm];
			h = rv1[k];
			f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g = vcn_math_hypo(f, 1.0);
			f = ((x-z)*(x+z) + 
			     h*((y/(f + (f >= 0.0 ? fabs(g):-fabs(g))))-h))/x;
			c = s = 1.0;
			for(j=l; j<=nm; j++){
				i = j+1;
				g = rv1[i];
				y = w[i];
				h = s*g;
				g = c*g;
				z = vcn_math_hypo(f, h);
				rv1[j] = z;
				c = f/z;
				s = h/z;
				f = x*c+g*s;
				g = g*c-x*s;
				h = y*s;
				y *= c;
				for(jj=0; jj<n; jj++){
					x = V[jj*n+j];
					z = V[jj*n+i];
					V[jj*n+j] = x*c+z*s;
					V[jj*n+i] = z*c-x*s;
				}
				z = vcn_math_hypo(f, h);
				w[j] = z;
				if(z){
					z = 1.0/z;
					c = f*z;
					s = h*z;
				}
				f = c*g+s*y;
				x = c*y-s*g;
				for(jj=0; jj<m; jj++){
					y = A[jj*n+j];
					z = A[jj*n+i];
					A[jj*n+j] = y*c+z*s;
					A[jj*n+i] = z*c-y*s;
				}
			}
			rv1[l] = 0.0;
			rv1[k] = f;
			w[k] = x;
		}
	}
	/* Free memory */
	free(rv1);
}


void vcn_matrix_svd_solve(const double *const U,
			  const double *const w, 
			  const double *const V, 
			  double *x,             /* Out */
			  const double *const b, 
			  int n, int m){
	/* Numerical Recipes in C.
	 * Solves Ax = b for a vector X, where A is specified by the arrays u, w, v 
	 * as returned by svd_dcmp, m and n are the dimensions of a, and will be equal
	 * for square matrices. b is the input right-hand side. x is the output solution 
	 * vector. No input quantities are destroyed, so the routine may be called 
	 * sequentially with different b's. 
	 */
	int jj, j, i;
	double s, *tmp;
	tmp = (double*)malloc(n*sizeof(double));
	for(j=0; j<n; j++){
		s = 0.0;
		if(w[j]){
			for(i=0; i<m; i++) s += U[i*n+j]*b[i];
			s /= w[j];
		}
		tmp[j] = s;
	}
	for(j=0; j<n; j++){
		s = 0.0;
		for(jj=0; jj<n; jj++) s += V[j*n+jj]*tmp[jj];
		x[j] = s;
	}
	free(tmp);
}


void vcn_matrix_forward_solve(const double *const L,
			      const double *const b, 
			      double *_x, uint32_t N)
/* Solve the system Lx = b, where L is a lower triangular matrix */
{
	/* This solver cannot be parallelized to guaranty the solution */
	for(uint32_t i=0; i< N; i++){
		_x[i] = b[i];
		for(uint32_t k=0; k < i; k++)
			_x[i] -= L[i*N+k]*_x[k];
		_x[i] /= L[i*N+i];
	}
}

void vcn_matrix_backward_solve(const double *const U,
			       const double *const b,
			       double *_x, uint32_t N)
/* Solve the system Ux = b, where U is a upper triangular matrix */
{
	/* This solver cannot be parallelized to guaranty the solution */
	for(int i=N-1; i>=0 ; i--){
		_x[i] = b[i];
		for(int k=i+1; k<N; k++)
			_x[i] -= U[i*N+k]*_x[k];
		_x[i] /= U[i*N+i];
	}
}

void vcn_sparse_eigen_power(const vcn_sparse_t* const A, int h,
			    double **_eigenvecs,/* Out */ 
			    double *_eigenvals, /* Out */
			    int *it,            /* Out */
			    double tolerance,
			    uint32_t omp_parallel_threads){
	/* The program must receive all the pointers allocated, where
	 *  > A is a vcn_sparse_t matrix
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
	double *p = (double*)calloc(A->N, sizeof(double));

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
	free(p);
}

int vcn_sparse_eigen_ipower(const vcn_sparse_t *const A,
			    int solver_type,
			    int h, double mu, 
			    double **_eigenvecs,/* Out */ 
			    double *_eigenvals, /* Out */
			    int* it,            /* Out */
			    double tolerance,
			    uint32_t omp_parallel_threads)
{
	/* The program must receive all the pointers allocated, where
	 *  > A is a vcn_sparse_t matrix
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
	double* p = (double*) calloc(A->N, sizeof(double));
	double* z = (double*) calloc(A->N, sizeof(double));

	/* Set M in A (copy ptr to modify const A, it will be restored) */
	vcn_sparse_t* A_ptr_copy = (vcn_sparse_t*)A;
	if (mu != 0.0)
		for (c = 0; c < A->N; c++)
			vcn_sparse_add(A_ptr_copy, c, c, -mu);            /* M = A - mu*I */

	/* LU Decomposition in case of LU Solver */
	vcn_sparse_t *L = NULL;
	vcn_sparse_t *U = NULL;
	if (solver_type == NB_SOLVER_CHK) {
		vcn_sparse_alloc_LU(A, &L, &U);
		vcn_sparse_decompose_Cholesky(A, L, U, omp_parallel_threads);
	} else if (solver_type == NB_SOLVER_LUD) {
		vcn_sparse_alloc_LU(A, &L, &U);
		vcn_sparse_decompose_LU(A, L, U, omp_parallel_threads);
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
			if (solver_type == NB_SOLVER_CHK || solver_type == NB_SOLVER_LUD)
				vcn_sparse_solve_LU(L, U, _eigenvecs[i], p);
			else if (solver_type == NB_SOLVER_CGJ)
				vcn_sparse_solve_CG_precond_Jacobi(A,_eigenvecs[i], p, 
								   vcn_sparse_get_size(A)*10, 1e-3, 
								   NULL, NULL,
								   omp_parallel_threads);
			else
				vcn_sparse_solve_conjugate_gradient(A,_eigenvecs[i], p, 
								    vcn_sparse_get_size(A)*10, 1e-3, 
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
			vcn_sparse_add(A_ptr_copy, c, c, mu);  /* A = M + mu*I */
    
	/* Destroy LU decomposition */
	if (solver_type == NB_SOLVER_CHK || solver_type == NB_SOLVER_LUD) {
		vcn_sparse_destroy(U);
		vcn_sparse_destroy(L);
	}
	/* Free memory */
	free(p);
	free(z);

	return 0;
}


void vcn_sparse_eigen_lanczos(const vcn_sparse_t* const A,
			      double *_eigenmax,/* Out */ 
			      double *_eigenmin,/* Out */
			      int* it,          /* Out */
			      double tolerance,
			      uint32_t omp_parallel_threads)
{
	/* The program must receive all the pointers allocated, where
	 *  > A is a vcn_sparse_t matrix
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
	double* alpha = (double*) calloc(A->N, sizeof(double));
	double* beta = (double*) calloc(A->N, sizeof(double));
	double* v = (double*) calloc(A->N, sizeof(double));
	double* w = (double*) calloc(A->N, sizeof(double));
	double* v_prev = (double*) calloc(A->N, sizeof(double));

	*it = 0;
	beta[*it] = 0;
	v[0] = 1;
	/* Start iterations */
	double normw2 = 1;
	double delta_eigen = 1;

	while(delta_eigen >= POW2(tolerance) &&
	      normw2 >= tolerance && *it < A->N){
		/* Step 1 and 2 */
		double ak = 0;
#pragma omp parallel for reduction(+:ak) num_threads(omp_parallel_threads) private(i, j)
		for(i=0; i<A->N; i++){
			w[i] = 0;
			for(j=0; j<A->rows_size[i]; j++)
				w[i] += A->rows_values[i][j]*v[A->rows_index[i][j]];
			w[i] -= beta[*it]*v_prev[i];
			ak += w[i]*v[i];
		}
		alpha[*it] = ak;
		/* Step 3 and 4 */
		*it = *it+1;
		normw2 = 0;
		for(i=0; i < A->N; i++){
			w[i] -= ak*v[i];
			normw2 += w[i]*w[i];
		}
		normw2 = sqrt(normw2);
		beta[*it] = normw2;
		/* Step 5 */
		for(i=0; i < A->N; i++){
			v_prev[i] = v[i];
			v[i] = w[i]/beta[*it];
		}
		/* Step 6 and 7 */
		if(*it > 1){
			double delta_eig1 = *_eigenmax;
			double delta_eigk = *_eigenmin;
			vcn_sparse_eigen_givens(alpha, beta, 1, _eigenmax, tolerance, *it);
			vcn_sparse_eigen_givens(alpha, beta, *it, _eigenmin, tolerance, *it);
      
			delta_eigen = fabs(delta_eig1-*_eigenmax);
			double tmp = fabs(delta_eigk-*_eigenmin);
			if(tmp > delta_eigen)
				delta_eigen = tmp;
		}else{
			*_eigenmax = ak;
			*_eigenmin = ak;
		}
	}

	/* Free memory */
	free(alpha);
	free(beta);
	free(w);
	free(v);
	free(v_prev);
}

void vcn_sparse_eigen_givens(const double* const main_diag, 
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
	double* p = (double*)calloc(N+1, sizeof(double));

	/* Init algorithm */
	int k = N;
	double a = main_diag[0]-fabs(uplw_diag[1]);
	double tmp = main_diag[k-1]-fabs(uplw_diag[k-1]);
	if(tmp < a)
		a = tmp;
	double b = main_diag[0]+fabs(uplw_diag[1]);
	tmp = main_diag[k-1]+fabs(uplw_diag[k-1]);
	if(tmp > b)
		b = tmp;
	for(l=1; l<k-1; l++){
		tmp = main_diag[l]-
			fabs(uplw_diag[l+1])-fabs(uplw_diag[l]);
		if(tmp < a)
			a = tmp;
		tmp = main_diag[l]+
			fabs(uplw_diag[l+1])+fabs(uplw_diag[l]);
		if(tmp > b)
			b = tmp;
	}
	/* Init iterations */
	while(fabs(b-a) > (fabs(a)+fabs(b)) * tolerance){
		/* Step 1 */
		*_eigenvalue = (a+b)/2;
		/* Step 2 */
		double r = 0;
		double s = 0;
		p[0] = 1;
		p[1] = main_diag[0] - *_eigenvalue;
		for(l=1; l<k; l++){
			p[l+1] = (main_diag[l]-*_eigenvalue)*p[l]-
				(uplw_diag[l])*(uplw_diag[l])*p[l-1];
		}
		for(l=1; l<k+1; l++){
			if(p[l]*p[l-1] <= 0)
				r++;
			if(p[l] == 0)
				s++;
		}
		double gamma = r-s;
		/* Step 3 */
		if(gamma > k-i)
			b = *_eigenvalue;
		else
			a = *_eigenvalue;
	}

	/* Free memory */
	free(p);
}

/* Meta-functions */
int meta_compare_data_bycol_2(const void *a, const void *b){
	/* The "2" anexed to the name is because the function is 
	 * already defined in the vcn_sparse_t struct implementation
	 */
	if((*(double**)a)[0] == (*(double**)b)[0])
		return (*(double**)a)[1] - (*(double**)b)[1];
	else
		return (*(double**)a)[0] - (*(double**)b)[0];

}
void vcn_sparse_read_mat4(vcn_sparse_t *A, const char *url, char *label){
	/* Read a vcn_sparse_t matrix named [label] from matlab v4 file */
 
	/* Open file */
	FILE *pfile = fopen(url,"rb");
	if(pfile == NULL){
		printf("ERROR: Impossible to open matlab v4 file.\n");
		exit(1);
	}
	/* Init int32[5] to store info */
	int32_t info[5];
	short found = 0;
	while(fread(info, 4, 5, pfile) != 0 && found != 1){
		uint32_t i;
		/* info[0] <- Type of structure and type of data
		 *  __________________________
		 *  |   0 := Complete matrix |
		 *  |   1 := Text            |
		 *  |   2 := Sparse matrix   |
		 *  |  +0 := double          |
		 *  | +10 := float           |
		 *  | +20 := int             |
		 *  | +30 := short           |
		 *  | +40 := unsigned short  |
		 *  | +50 := unsigned char   |
		 *  |________________________|
		 *
		 * info[4] <- Name length + 1
		 */
		char* name = malloc(info[4] * sizeof(char));
		/* Read name */
		if(fread(name, sizeof(char), info[4], pfile) == 0){
			printf("ERROR: The matlab v4 file is corrupted.\n");
			exit(1);
		}
		if(info[0] == 0){
			/* Jump to the next structure */
			if(fseek(pfile, info[1]*info[2]*sizeof(double), SEEK_CUR) != 0){
				printf("ERROR: The matlab v4 file is corrupted.\n");
				exit(1);
			}
		}else if(info[0] == 2){
			/* Vcn_Sparse_T matrix
			 *
			 *  info[1] <- Number of non-zero entries (nnz) +1
			 *  info[2] <- Type of values {0(Real) | 1(Complex)}
			 *  info[3] <- Zero
			 */
			if(strcmp(name,label) == 1){
				/* If "name" and "label" are not equal
				 * then jump to the next structure 
				 */
				if(fseek(pfile, (3*(info[1]-1)+3)*sizeof(double), SEEK_CUR) != 0){
					printf("ERROR: The matlab v4 file is corrupted.\n");
					exit(1);
				}
			}else{
				found = 1;

				uint32_t nnz = (uint32_t)(info[1]-1);
				/* Verification of use of real numbers */
				if(info[2] == 4){
					printf("ERROR: Complex numbers unsupported.\n");
					exit(1);
				}
				double *irows = (double*)malloc(nnz*sizeof(double));
				double *icols = (double*)malloc(nnz*sizeof(double));
				double *values = (double*)malloc(nnz*sizeof(double));
				double N;
				/* Read row's index */
				if(fread(irows, sizeof(double), nnz, pfile) == 0){
					printf("ERROR: The matlab v4 file is corrupted.\n");
					exit(1);
				}
				/* Read number of rows */
				if(fread(&N, sizeof(double), 1, pfile) == 0){
					printf("ERROR: The matlab v4 file is corrupted.\n");
					exit(1);
				}
				A->N = (uint32_t)N;
				A->rows_values = (double**)malloc(A->N*sizeof(void*));
				A->rows_index = (uint32_t**)malloc(A->N*sizeof(void*));
				A->rows_size = (uint32_t*)calloc(A->N,sizeof(uint32_t));

				/* Read col's index */
				if(fread(icols, sizeof(double), nnz, pfile) == 0){
					printf("ERROR: The matlab v4 file is corrupted.\n");
					exit(1);
				}
				/* Read number of rows */
				if(fread(&N, sizeof(double), 1, pfile) == 0){
					printf("ERROR: The matlab v4 file is corrupted.\n");
					exit(1);
				}
				A->N = (uint32_t)N;

				/* Read values */
				if(fread(values, sizeof(double), nnz, pfile) == 0){
					printf("ERROR: The matlab v4 file is corrupted.\n");
					exit(1);
				}

				uint32_t* rows_icol = (uint32_t*)calloc(N,sizeof(uint32_t));
				for(i=0; i< nnz; i++)
					A->rows_size[(uint32_t)irows[i]-1]++;
	
				for(i=0; i< N; i++){
					A->rows_index[i] = (uint32_t*)calloc(A->rows_size[i],sizeof(uint32_t));
					A->rows_values[i] = (double*)calloc(A->rows_size[i],sizeof(double));
				}

				for(i=0; i<nnz; i++){
					uint32_t irow = irows[i]-1;
					A->rows_index[irow][rows_icol[irow]] = (uint32_t)icols[i]-1;
					A->rows_values[irow][rows_icol[irow]] = values[i];
					rows_icol[irow] ++;
				}

				/* Sort data by columns */
				for(i=0; i<A->N; i++)
					vcn_qsort(A->rows_index[i], A->rows_size[i], 
						  sizeof(uint32_t), vcn_compare_uint32);

				/* Free memory */
				free(irows);
				free(icols);
				free(values);
				free(rows_icol);
			}
		}else{
			printf("ERROR: Support only for complete and sparse matrix with double precision.\n");
			exit(1);
		}
		/* Free memory */
		free(name);
	}
	/* Close file */
	fclose(pfile);

	if(found != 1){
		printf("ERROR: Sparse matrix \"%s\" not found in \"%s\".\n",label,url);
		exit(1);
	}
}

void vcn_sparse_save_mat4(const vcn_sparse_t *const A,
			  const char *url, char *label)
/* Write a sparse matrix named [label] in matlab v4 format */
{
	uint32_t i, j;
	/* Open file */
	FILE *pfile = fopen(url,"ab");
	if(pfile == NULL){
		printf("ERROR: Impossible to open matlab v4 file.\n");
		exit(1);
	}
	uint32_t nnz = vcn_sparse_get_nnz(A);
	int32_t info[5];
	info[0] = 2; /* Sparse matrix with double precision */
	info[1] = nnz + 1;
	info[2] = 3; /* Real values */
	info[3] = 0; /* Strictly zero by matlab v4 definition*/
	int name_length = strlen(label);
	info[4] = name_length+1;
	/* Write headers of struct */
	fwrite(info, 4, 5, pfile);
	/* Write label of struct */
	fwrite(label, sizeof(char), name_length, pfile);
	char end_label = '\0';
	fwrite(&end_label, sizeof(char), 1, pfile);
	/* Write sparse matrix data */
	double** data = (double**)malloc(nnz*sizeof(void*));
	for(i=0; i<nnz; i++)
		data[i] = (double*)malloc(3*sizeof(double));
	uint32_t idata = 0;
	for(i=0; i< A->N; i++){
		for(j=0; j< A->rows_size[i]; j++){
			data[idata][0] = A->rows_index[i][j];
			data[idata][1] = i;
			data[idata][2] = A->rows_values[i][j];
			idata++;
		}
	}
	/* Sort vectors by column and then by rows for octave compability */
	qsort(data, nnz, sizeof(double*), meta_compare_data_bycol_2);/* TEMPORAL: Use nb_qsort */
	/* Convert data to vectors to write in octave file */
	double* irows = (double*)malloc(nnz*sizeof(double));
	double* icols = (double*)malloc(nnz*sizeof(double));
	double* values = (double*)malloc(nnz*sizeof(double));
	for(i=0; i < nnz; i++){
		icols[i] = data[i][0]+1;
		irows[i] = data[i][1]+1;
		values[i] = data[i][2];
	}
	/* Write row's index */
	fwrite(irows, sizeof(double), nnz, pfile);
	/* Write number of rows */
	double N = (double)A->N;
	fwrite(&N, sizeof(double), 1, pfile);
	/* Write col's index */
	fwrite(icols, sizeof(double), nnz, pfile);
	/* Write number of cols */
	fwrite(&N, sizeof(double), 1, pfile);
	/* Write values */
	fwrite(values, sizeof(double), nnz, pfile);
	/* Write a zero value for octave compability */
	double zero = 0;
	fwrite(&zero, sizeof(double), 1, pfile);
	/* Close file */
	fclose(pfile);
}

/* Functions to load and write Octave files (MAT-File 4 format) */
void vcn_mat4_printf(const char *url)
{
	/* Read the octave file to get and print information 
	 * about stored objects.
	 */
	/* Open file */
	FILE *pfile = fopen(url,"rb");
	if(pfile == NULL){
		printf("ERROR: Impossible to open matlab v4 file.\n");
		exit(1);
	}
	/* Init int32[5] to store info */
	int32_t info[5];
	/* Read items */
	while(fread(info, 4, 5, pfile) != 0){
		/* info[0] <- Type of structure and type of data
		 *  __________________________
		 *  |   0 := Complete matrix |
		 *  |   1 := Text            |
		 *  |   2 := Vcn_Sparse_T matrix   |
		 *  |  +0 := double          |
		 *  | +10 := float           |
		 *  | +20 := int             |
		 *  | +30 := short           |
		 *  | +40 := unsigned short  |
		 *  | +50 := unsigned char   |
		 *  |________________________|
		 *
		 * info[4] <- Name length + 1
		 */
		char* name = malloc(info[4]*sizeof(char));
		/* Read name */
		if(fread(name, sizeof(char), info[4], pfile) == 0){
			printf("ERROR: The matlab v4 file is corrupted.\n");
			exit(1);
		}
		if(info[0] == 0){
			/* Complete matrix (and vectors seen as matrix by octave)
			 *
			 *  info[1] <- Number of rows
			 *  info[2] <- Number of cols
			 *  info[3] <- Type of values {0(Real) | 1(Complex)}
			 */

			/* Verification of use of real numbers */
			if(info[3] == 1){
				printf("ERROR: Complex numbers unsupported.\n");
				exit(1);
			}
			/* Show information about structure */
			if(info[1] == 1)
				printf("%s <- Vector (size: %d)\n", name, info[2]);
			else if(info[2] == 1)
				printf("%s <- Vector (size: %d)\n", name, info[1]);
			else
				printf("%s <- Complete matrix (size: %dx%d)\n", name, info[1],info[2]);
			/* Jump to the next structure */
			if(fseek(pfile, info[1]*info[2]*sizeof(double), SEEK_CUR) != 0){
				printf("ERROR: The matlab v4 file is corrupted.\n");
				exit(1);
			}
		}else if(info[0] == 2){
			/* Vcn_Sparse_T matrix
			 *
			 *  info[1] <- Number of non-zero entries (nnz) +1
			 *  info[2] <- Type of values {0(Real) | 1(Complex)}
			 *  info[3] <- Zero
			 */
			uint32_t nnz = info[1]-1;
			/* Verification of use of real numbers */
			if(info[2] == 1){
				printf("ERROR: Complex numbers unsupported.\n");
				exit(1);
			}
			double N;
			/* Jump to get number of rows */
			if(fseek(pfile, nnz*sizeof(double), SEEK_CUR) != 0){
				printf("ERROR: The matlab v4 file is corrupted.\n");
				exit(1);
			}
			/* Read number of rows */
			if(fread(&N, sizeof(double), 1, pfile) == 0){
				printf("ERROR: The matlab v4 file is corrupted.\n");
				exit(1);
			}
			/* Jump to get number of cols */
			if(fseek(pfile, nnz*sizeof(double), SEEK_CUR) != 0){
				printf("ERROR: The matlab v4 file is corrupted.\n");
				exit(1);
			}
			/* Read number of rows */
			if(fread(&N, sizeof(double), 1, pfile) == 0){
				printf("ERROR: The matlab v4 file is corrupted.\n");
				exit(1);
			}
			/* Show information about structure */
			printf("%s <- Vcn_Sparse_T matrix (size: %ix%i, nnz: %d).\n",
			       name, (uint32_t)N, (uint32_t)N, info[1]);
			/* Jump to the next structure */
			if(fseek(pfile, (nnz+1)*sizeof(double), SEEK_CUR) != 0){
				printf("ERROR: The matlab v4 file is corrupted.\n");
				exit(1);
			}
		}else{
			printf("ERROR: Support only for complete and vcn_sparse_t matrix with double precision.\n");
			exit(1);
		}
		/* Free memory */
		free(name);
	}
	/* Close file */
	fclose(pfile);
}

short vcn_mat4_exist(const char *url, char* label)
{
	/* Search a structure in octave file.
	 * Return 1 if found the structure and 0 if not.
	 */
	short found = 0;
	/* Open file */
	FILE *pfile = fopen(url,"rb");
	if(pfile == NULL){
		printf("ERROR: Impossible to open matlab v4 file.\n");
		exit(1);
	}
	/* Init int32[5] to store info */
	int32_t info[5];
	/* Read items */
	while(fread(info, 4, 5, pfile) != 0 && found == 0){
		/* info[0] <- Type of structure and type of data
		 *  __________________________
		 *  |   0 := Complete matrix |
		 *  |   1 := Text            |
		 *  |   2 := Vcn_Sparse_T matrix   |
		 *  |  +0 := double          |
		 *  | +10 := float           |
		 *  | +20 := int             |
		 *  | +30 := short           |
		 *  | +40 := unsigned short  |
		 *  | +50 := unsigned char   |
		 *  |________________________|
		 *
		 * info[4] <- Name length + 1
		 */
		char* name = malloc(info[4]*sizeof(char));
		/* Read name */
		if(fread(name, sizeof(char), info[4], pfile) == 0){
			printf("ERROR: The matlab v4 file is corrupted.\n");
			exit(1);
		}
		if(strcmp(name,label) == 0)
			found = 1;
		else if(info[0] == 0){
			/* Jump to the next structure */
			if(fseek(pfile, info[1]*info[2]*sizeof(double), SEEK_CUR) != 0){
				printf("ERROR: The matlab v4 file is corrupted.\n");
				exit(1);
			}
		}else if(info[0] == 2){
			/* Jump to the next structure */
			if(fseek(pfile, (3*(info[1]-1)+3)*sizeof(double), SEEK_CUR) != 0){
				printf("ERROR: The matlab v4 file is corrupted.\n");
				exit(1);
			}
		}else{
			printf("ERROR: Support only for complete and vcn_sparse_t matrix with double precision.\n");
			exit(1);
		}
		/* Free memory */
		free(name);
	}
	/* Close file */
	fclose(pfile);

	return found;
}

void vcn_mat4_clear(const char *url)
{
	FILE *pfile = fopen(url,"wb");
	fclose(pfile);
}

void vcn_mat4_read_vec(const char *url, char *label, double *_x)
{
	/* Read a vector named [label] from matlab v4 file */
 
	/* Open file */
	FILE *pfile = fopen(url,"rb");
	if(pfile == NULL){
		printf("ERROR: Impossible to open matlab v4 file.\n");
		exit(1);
	}
	/* Init int32[5] to store info */
	int32_t info[5];
	short found = 0;
	while(fread(info, 4, 5, pfile) != 0 && found != 1){
		/* info[0] <- Type of structure and type of data
		 *  __________________________
		 *  |   0 := Complete matrix |
		 *  |   1 := Text            |
		 *  |   2 := Sparse matrix   |
		 *  |  +0 := double          |
		 *  | +10 := float           |
		 *  | +20 := int             |
		 *  | +30 := short           |
		 *  | +40 := unsigned short  |
		 *  | +50 := unsigned char   |
		 *  |________________________|
		 *
		 * info[4] <- Name length + 1
		 */
		char* name = malloc(info[4]*sizeof(char));
		/* Read name */
		if(fread(name, sizeof(char), info[4], pfile) == 0){
			printf("ERROR: The matlab v4 file is corrupted.\n");
			exit(1);
		}
		if(info[0] == 0){
			/* Vector seen as matrix by octave
			 *
			 *  info[1] <- Number of rows
			 *  info[2] <- Number of cols
			 *  info[3] <- Type of values {0(Real) | 1(Complex)}
			 */
			int min_size = MIN(info[1], info[2]);
			if(strcmp(name,label) == 1 || min_size > 1){
				/* If "name" and "label" are not equal
				 * then jump to the next structure 
				 */
				if(fseek(pfile, info[1]*info[2]*sizeof(double), SEEK_CUR) != 0){
					printf("ERROR: The matlab v4 file is corrupted.\n");
					exit(1);
				}
			}else{
				/* Verification of use of real numbers */
				if(info[3] == 1){
					printf("ERROR: Complex numbers unsupported.\n");
					exit(1);
				}
				found = 1;
				/* Initialize vector */
				uint32_t size;
				if(info[1] == 1)
					size = info[2];
				else
					size = info[1];
				_x = (double*)calloc(size, sizeof(double));
				if(fread(_x, sizeof(double), size, pfile) == 0){
					printf("ERROR: The matlab v4 file is corrupted.\n");
					exit(1);
				}
			}
		}else if(info[0] == 2){
			/* Jump to the next structure */
			if(fseek(pfile, (3*(info[1]-1)+3)*sizeof(double), SEEK_CUR) != 0){
				printf("ERROR: The matlab v4 file is corrupted.\n");
				exit(1);
			}
		}else{
			printf("ERROR: Support only for complete and sparse matrix with double precision.\n");
			exit(1);
		}
		/* Free memory */
		free(name);
	}
	/* Close file */
	fclose(pfile);

	if(found != 1){
		printf("ERROR: Vector \"%s\" not found in \"%s\".\n",label,url);
		exit(1);
	}
}

void vcn_mat4_save_vec(const char *url, char *label, 
		       const double *const x, uint32_t N){
	/* Write a sparse matrix named [label] in matlab v4 format */

	/* Open file */
	FILE *pfile = fopen(url,"ab");
	if(pfile == NULL){
		printf("ERROR: Impossible to open matlab v4 file.\n");
		exit(1);
	}
	int32_t info[5];
	info[0] = 0; /* Complete matrix with double precision */
	info[1] = N;     /* N */
	info[2] = 1;     /* N */
	info[3] = 0; /* Real values */
	int name_length = strlen(label);
	info[4] = name_length+1;
	/* Write headers of struct */
	fwrite(info, 4, 5, pfile);
	/* Write label of struct */
	fwrite(label, sizeof(char), name_length, pfile);
	char end_label = '\0';
	fwrite(&end_label, sizeof(char), 1, pfile);
	/* Write vector values */
	fwrite(x, sizeof(double), N, pfile);
	/* Close file */
	fclose(pfile);
}

void vcn_mat4_read_mtx(const char *url, char *label, double *_A){
	/* Read a vector named [label] from matlab v4 file */
 
	/* Open file */
	FILE *pfile = fopen(url,"rb");
	if(pfile == NULL){
		printf("ERROR: Impossible to open matlab v4 file.\n");
		exit(1);
	}
	/* Init int32[5] to store info */
	int32_t info[5];
	short found = 0;
	while(fread(info, 4, 5, pfile) != 0 && found != 1){
		/* info[0] <- Type of structure and type of data
		 *  __________________________
		 *  |   0 := Complete matrix |
		 *  |   1 := Text            |
		 *  |   2 := Sparse matrix   |
		 *  |  +0 := double          |
		 *  | +10 := float           |
		 *  | +20 := int             |
		 *  | +30 := short           |
		 *  | +40 := unsigned short  |
		 *  | +50 := unsigned char   |
		 *  |________________________|
		 *
		 * info[4] <- Name length + 1
		 */
		char* name = malloc(info[4]*sizeof(char));

		/* Read name */
		if(fread(name, sizeof(char), info[4], pfile) == 0){
			printf("ERROR: The matlab v4 file is corrupted.\n");
			exit(1);
		}

		if(info[0] == 0){
			/* Vector seen as matrix by octave
			 *
			 *  info[1] <- Number of rows
			 *  info[2] <- Number of cols
			 *  info[3] <- Type of values {0(Real) | 1(Complex)}
			 */
			if(strcmp(name,label) == 1){
				/* If "name" and "label" are not equal
				 * then jump to the next structure 
				 */
				if(fseek(pfile, info[1]*info[2]*sizeof(double), SEEK_CUR) != 0){
					printf("ERROR: The matlab v4 file is corrupted.\n");
					exit(1);
				}
			}else{
				/* Verification of use of real numbers */
				if(info[3] == 1){
					printf("ERROR: Complex numbers unsupported.\n");
					exit(1);
				}
				found = 1;
				/* Initialize vector */
				uint32_t N = info[1];
				_A = (double*)calloc(N*N,sizeof(double));
				if(fread(_A, sizeof(double), N*N, pfile) == 0){
					printf("ERROR: The matlab v4 file is corrupted.\n");
					exit(1);
				}
			}
		}else if(info[0] == 2){
			/* Jump to the next structure */
			if(fseek(pfile, (3*(info[1]-1)+3)*sizeof(double), SEEK_CUR) != 0){
				printf("ERROR: The matlab v4 file is corrupted.\n");
				exit(1);
			}
		}else{
			printf("ERROR: Support only for complete and sparse matrix with double precision.\n");
			exit(1);
		}
		/* Free memory */
		free(name);
	}
	/* Close file */
	fclose(pfile);

	if(found != 1){
		printf("ERROR: Complete matrix \"%s\" not found in \"%s\".\n",label,url);
		exit(1);
	}
}

void vcn_mat4_write_mtx(const char *url, char *label,
			const double *const A, uint32_t N){
	/* Write a matrix named [label] in matlab v4 format */

	/* Open file */
	FILE *pfile = fopen(url,"ab");
	if(pfile == NULL){
		printf("ERROR: Impossible to open matlab v4 file.\n");
		exit(1);
	}
	int32_t info[5];
	info[0] = 0; /* Complete matrix with double precision */
	info[1] = N;
	info[2] = N;
	info[3] = 0; /* Real values */
	int name_length = strlen(label);
	info[4] = name_length+1;
	/* Write headers of struct */
	fwrite(info, 4, 5, pfile);
	/* Write label of struct */
	fwrite(label, sizeof(char), name_length, pfile);
	char end_label = '\0';
	fwrite(&end_label, sizeof(char), 1, pfile);
	/* Transpose matrix */
	double* At = (double*) calloc(N*N, sizeof(double));
	uint32_t i, j;
	for(i=0; i < N; i++)
		for(j=0; j < N; j++)
			At[i*N + j] = A[j*N + i];

	/* Write matrix values */
	fwrite(At, sizeof(double), N*N, pfile);
	free(At);
	/* Close file */
	fclose(pfile);
}
