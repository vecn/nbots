#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "nb/math_bot.h"
#include "nb/memory_bot.h"
#include "nb/container_bot/array.h"
#include "nb/solver_bot/sparse.h"

#include "sparse_struct.h"

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

static inline double get_norm(double* x, uint32_t N)
{
	double n = 0;
	for(uint32_t i=0; i<N; i++)
		n += x[i]*x[i];
	return sqrt(n);
}

static inline vcn_sparse_t* sparse_allocate(uint32_t N)
{
	vcn_sparse_t* A = nb_allocate_mem(sizeof(*A));
	A->rows_values = nb_allocate_mem(N * sizeof(*(A->rows_values)));
	A->rows_index = nb_allocate_mem(N * sizeof(*(A->rows_index)));
	A->rows_size = nb_allocate_zero_mem(N * sizeof(*(A->rows_size)));
	A->N = N;
	return A;
}

vcn_sparse_t* vcn_sparse_create(const nb_graph_t *const restrict graph,
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
		nb_free_mem(A->rows_values[i]);
		nb_free_mem(A->rows_index[i]);
	}
	nb_free_mem(A->rows_values);
	nb_free_mem(A->rows_index);
	nb_free_mem(A->rows_size);
	nb_free_mem(A);
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
	double* br = (double*)nb_allocate_mem(N * sizeof(double));
	for(uint32_t i=0; i < N; i++)
		br[i] = b[perm[i]];
	return br;
}

void vcn_sparse_get_transpose(const vcn_sparse_t *A, vcn_sparse_t *_At)
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

void vcn_sparse_transpose(vcn_sparse_t *A)
{
	for (uint32_t i = 0; i < A->N; i++) {
		for (uint32_t q = 0; A->rows_index[i][q] < i; q++) {
			uint32_t j = A->rows_index[i][q];
			uint32_t jc = sparse_bsearch_row(A, j, i, 0,
							 A->rows_size[j]-1);
			double aux = A->rows_values[i][q];
			A->rows_values[i][q] = A->rows_values[j][jc];
			A->rows_values[j][jc] = aux;
		}
	}
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

double vcn_sparse_get_usym(const vcn_sparse_t *const A)
{
	double sum = 0;
	for (uint32_t i = 0; i < A->N; i++) {
		for (uint32_t j = 0; A->rows_index[i][j] < i; j++) {
			double Aij = A->rows_values[i][j];
			uint32_t k = A->rows_index[i][j];
			uint32_t l = sparse_bsearch_row(A, k, i, 0,
							A->rows_size[k]-1);
			double Aji = A->rows_values[k][l];
			sum += POW2(Aij - Aji); 
		}
	}
	return sqrt(sum);
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

void vcn_sparse_multiply_vector(const vcn_sparse_t* A, const double* in,
				double* out, uint32_t omp_parallel_threads)
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
//	float* cells = nb_allocate_zero_mem(POW2(plot_size) * sizeof(float));
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
//	row_pointers = png_nb_allocate_mem(png_ptr, h*sizeof(png_byte*));
//	for(i=0; i<h; i++){
//		png_byte *row=
//			png_nb_allocate_mem(png_ptr, sizeof(uint8_t)*w*pixel_size);
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
//		png_nb_free_mem(png_ptr, row_pointers[i]);
//	png_nb_free_mem(png_ptr, row_pointers);
//	png_destroy_write_struct(&png_ptr, &info_ptr);
//
//	nb_free_mem(cells);
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
	double* c = (double*)nb_allocate_mem(vcn_sparse_get_size(A)*sizeof(double));
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
	nb_free_mem(c);

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
	char *memblock = nb_allocate_zero_mem(5 * vec_size);	
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
		for (uint32_t i = 0; i < A->N; i++) {
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
		for (uint32_t i=0; i< A->N; i++) {
			_x[i] += alphak*p[i];
			g[i] += alphak*w[i];
			q[i] = g[i]/Aii[i];
			dot_gkqk += g[i]*q[i];
		}

		double betak = dot_gkqk/dot_gq;
		
#pragma omp parallel for num_threads(omp_parallel_threads)
		for (uint32_t i=0; i< A->N; i++)
			p[i] = -q[i]+betak*p[i];
		k++;
	}
	/* Free memory */
	nb_free_mem(memblock);

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
	double* g = nb_allocate_zero_mem(A->N * sizeof(double));
	double* p = nb_allocate_zero_mem(A->N * sizeof(double));
	double* q = nb_allocate_zero_mem(A->N * sizeof(double));
	double* w = nb_allocate_zero_mem(A->N * sizeof(double));

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
	nb_free_mem(g);
	nb_free_mem(p);
	nb_free_mem(q);
	nb_free_mem(w);
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
	double *D = nb_allocate_zero_mem(A->N * sizeof(double));
	double *siD = nb_allocate_zero_mem(A->N * sizeof(double));

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
		G->rows_index[i] = nb_allocate_zero_mem(isize *
							sizeof(uint32_t));
		G->rows_values[i] = nb_allocate_zero_mem(isize *
							 sizeof(double));

		Gt->rows_size[i] = isizet;
		Gt->rows_index[i] = nb_allocate_zero_mem(isizet *
							 sizeof(uint32_t));
		Gt->rows_values[i] = nb_allocate_zero_mem(isizet *
							  sizeof(double));
	}

#pragma omp parallel for num_threads(omp_parallel_threads)
	for(uint32_t i=0; i < A->N; i++){
		/* Compute values of ~G */
		double* subA =
			nb_allocate_zero_mem(POW2(G->rows_size[i]) *
					     sizeof(double));
		/* The data of vector g is not allocated, is a pointer to each row of ~G */
		double* subg = G->rows_values[i];
		double *delta = nb_allocate_zero_mem(G->rows_size[i] *
						     sizeof(double));
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
		double* L = nb_allocate_zero_mem(POW2(G->rows_size[i]) *
						 sizeof(double));
		vcn_matrix_cholesky_decomposition(subA, L, G->rows_size[i]);
		vcn_matrix_cholesky_solve(L, delta, subg, G->rows_size[i]);
		/* Finally do G = [~G]*D   */
		for(uint32_t q=0; q < G->rows_size[i]; q++)
			G->rows_values[i][q] *= D[G->rows_index[i][q]];

		/* Free memory */
		nb_free_mem(subA);
		nb_free_mem(L);
		nb_free_mem(delta);
	}
	/* Store G transposed */
	vcn_sparse_get_transpose(G,Gt);

	/* Free memory */
	nb_free_mem(D);
	nb_free_mem(siD);

	/* Solve Ax = b with Conjugate Gradient method */
	double* r = nb_allocate_zero_mem(A->N * sizeof(double));
	double* p = nb_allocate_zero_mem(A->N * sizeof(double));
	double* w = nb_allocate_zero_mem(A->N * sizeof(double));
	double* Gr = nb_allocate_zero_mem(A->N * sizeof(double));
	double* Mr = nb_allocate_zero_mem(A->N * sizeof(double));

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
	nb_free_mem(r);
	nb_free_mem(p);
	nb_free_mem(w);
	nb_free_mem(Gr);
	nb_free_mem(Mr);

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
	uint32_t *L_size = nb_allocate_zero_mem(A->N * sizeof(*L_size));

	sc_set** r = nb_allocate_zero_mem(A->N * sizeof(*r));

	for (uint32_t j = 0; j < A->N; j++) { 
		sc_set* lj = NULL;
		uint32_t lj_size = 1;
		uint32_t _i = sparse_bsearch_row(A, j, j, 0,
						 A->rows_size[j]-1);

		/* lj <- aj ************************************************/
		sc_set *iterator_lj;                                     /**/
		for(uint32_t i = _i+1; i<A->rows_size[j]; i++){     /**/
			sc_set* aji = (sc_set*)nb_allocate_mem(sizeof(sc_set));         /**/
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
							sc_set* node = (sc_set*)nb_allocate_mem(sizeof(sc_set)); /**/
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
					lj = (sc_set*)nb_allocate_mem(sizeof(sc_set));               /**/
					lj->index = _Lt->rows_index[i][k];                  /**/
					lj->next = NULL;                                    /**/
					lj_size++;                                          /**/
					k++;                                                /**/
				}                                                     /**/
				sc_set* iterator_lj = lj;                             /**/
				while(k < _Lt->rows_size[i]){                         /**/
					if(_Lt->rows_index[i][k] != j){                     /**/
						sc_set* node = (sc_set*)nb_allocate_mem(sizeof(sc_set));   /**/
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
						sc_set* node = (sc_set*)nb_allocate_mem(sizeof(sc_set));   /**/
						node->index = j;                                  /**/
						node->next = NULL;                                /**/
						iterator_rp->next = node;                         /**/
					}                                                   /**/
				}                                                     /**/
			}else{                                                  /**/
				r[p] = (sc_set*)nb_allocate_mem(sizeof(sc_set));               /**/
				r[p]->index = j;                                      /**/
				r[p]->next = NULL;                                    /**/
			}                                                       /**/
			/**********************************************************/
		}

		/*************** Allocate the jth row of "Lt" **********************/
		_Lt->rows_values[j] = nb_allocate_zero_mem(lj_size *
							   sizeof(double));   /**/
		_Lt->rows_index[j] = nb_allocate_zero_mem(lj_size *
							  sizeof(uint32_t));        /**/
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
			nb_free_mem(lj_free);                                                 /**/
		}                                                                /**/
		/*******************************************************************/
    
	}
	/****************************** Allocate "L" ****************************/
	for(uint32_t i=0; i<_L->N; i++){                                 /**/
		L_size[i]++;            /* To include main diagonal */              /**/
		_L->rows_values[i] = nb_allocate_zero_mem(L_size[i] *
							  sizeof(double));     /**/
		_L->rows_index[i] = nb_allocate_zero_mem(L_size[i] *
							 sizeof(uint32_t));          /**/
		_L->rows_size[i] = L_size[i];                                       /**/
	}                                                                     /**/
	/* Cycle to set columns in L */                                       /**/
	uint32_t* L_index = nb_allocate_zero_mem(_L->N * sizeof(uint32_t));                    /**/
	/* This "for" can not be parallelizad to guaranty increasing sort */  /**/
	for(uint32_t i=0; i<_Lt->N; i++){                                /**/
		for(uint32_t j=0; j<_Lt->rows_size[i]; j++){                   /**/
			uint32_t index = _Lt->rows_index[i][j];                               /**/
			_L->rows_index[index][L_index[index]++] = i;                      /**/
		}                                                                   /**/
	}                                                                     /**/
	/************************************************************************/

	/* Free memory */
	nb_free_mem(L_size);
	nb_free_mem(L_index);
	for(uint32_t i=0; i< A->N; i++){
		sc_set* iterator_ri = r[i];
		while(iterator_ri != NULL){
			sc_set* rm = iterator_ri;
			iterator_ri = (sc_set*)iterator_ri->next;
			nb_free_mem(rm);
		}
	}
	nb_free_mem(r);
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
	double* z = nb_allocate_zero_mem(L->N * sizeof(*z));
	vcn_sparse_forward_solve(L, b, z);
	vcn_sparse_backward_solve(U, z, _x);
	nb_free_mem(z);
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

int vcn_sparse_eigen_ipower(const vcn_sparse_t *const A,
			    nb_solver_t solver,
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
	double* p = nb_allocate_zero_mem(A->N * sizeof(double));
	double* z = nb_allocate_zero_mem(A->N * sizeof(double));

	/* Set M in A (copy ptr to modify const A, it will be restored) */
	vcn_sparse_t* A_ptr_copy = (vcn_sparse_t*)A;
	if (mu != 0.0)
		for (c = 0; c < A->N; c++)
			vcn_sparse_add(A_ptr_copy, c, c, -mu);            /* M = A - mu*I */

	/* LU Decomposition in case of LU Solver */
	vcn_sparse_t *L = NULL;
	vcn_sparse_t *U = NULL;
	if (NB_SOLVER_CHK == solver) {
		vcn_sparse_alloc_LU(A, &L, &U);
		vcn_sparse_decompose_Cholesky(A, L, U, omp_parallel_threads);
	} else if (NB_SOLVER_LUD == solver) {
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
			if (NB_SOLVER_CHK == solver ||
			    NB_SOLVER_LUD == solver)
				vcn_sparse_solve_LU(L, U, _eigenvecs[i], p);
			else if (NB_SOLVER_CGJ == solver)
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
	if (NB_SOLVER_CHK == solver || NB_SOLVER_LUD == solver) {
		vcn_sparse_destroy(U);
		vcn_sparse_destroy(L);
	}
	/* Free memory */
	nb_free_mem(p);
	nb_free_mem(z);

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
			vcn_sparse_eigen_givens(alpha, beta, 1, _eigenmax, tolerance, *it);
			vcn_sparse_eigen_givens(alpha, beta, *it, _eigenmin, tolerance, *it);
      
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
