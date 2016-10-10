#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "nb/math_bot.h"
#include "nb/memory_bot.h"
#include "nb/container_bot.h"

#include "nb/solver_bot/sparse/sparse.h"
#include "nb/solver_bot/matlab_v4.h"

#include "sparse/sparse_struct.h"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define POW2(a) ((a)*(a))

static int meta_compare_data_bycol(const void *a, const void *b);

void vcn_sparse_read_mat4(vcn_sparse_t *A, const char *url, char *label)
{
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
		char* name = nb_allocate_mem(info[4] * sizeof(char));
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
				double *irows = (double*)nb_allocate_mem(nnz*sizeof(double));
				double *icols = (double*)nb_allocate_mem(nnz*sizeof(double));
				double *values = (double*)nb_allocate_mem(nnz*sizeof(double));
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
				A->rows_values = (double**)nb_allocate_mem(A->N*sizeof(void*));
				A->rows_index = (uint32_t**)nb_allocate_mem(A->N*sizeof(void*));
				A->rows_size = (uint32_t*)nb_allocate_zero_mem(A->N * sizeof(uint32_t));

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

				uint32_t* rows_icol = nb_allocate_zero_mem(N * sizeof(uint32_t));
				for(i=0; i< nnz; i++)
					A->rows_size[(uint32_t)irows[i]-1]++;
	
				for(i=0; i< N; i++){
					A->rows_index[i] = nb_allocate_zero_mem(A->rows_size[i] * sizeof(uint32_t));
					A->rows_values[i] = nb_allocate_zero_mem(A->rows_size[i] * sizeof(double));
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
				nb_free_mem(irows);
				nb_free_mem(icols);
				nb_free_mem(values);
				nb_free_mem(rows_icol);
			}
		}else{
			printf("ERROR: Support only for complete and sparse matrix with double precision.\n");
			exit(1);
		}
		/* Free memory */
		nb_free_mem(name);
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
	double** data = (double**)nb_allocate_mem(nnz*sizeof(void*));
	for(i=0; i<nnz; i++)
		data[i] = (double*)nb_allocate_mem(3*sizeof(double));
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
	qsort(data, nnz, sizeof(double*), 
	      meta_compare_data_bycol);/* TEMPORAL: Use nb_qsort */
	/* Convert data to vectors to write in octave file */
	double* irows = nb_allocate_mem(nnz*sizeof(double));
	double* icols = nb_allocate_mem(nnz*sizeof(double));
	double* values = nb_allocate_mem(nnz*sizeof(double));
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

static int meta_compare_data_bycol(const void *a, const void *b)
{
	if((*(double**)a)[0] == (*(double**)b)[0])
		return (*(double**)a)[1] - (*(double**)b)[1];
	else
		return (*(double**)a)[0] - (*(double**)b)[0];
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
		char* name = nb_allocate_mem(info[4]*sizeof(char));
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
		nb_free_mem(name);
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
		char* name = nb_allocate_mem(info[4]*sizeof(char));
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
		nb_free_mem(name);
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
		char* name = nb_allocate_mem(info[4]*sizeof(char));
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
				_x = (double*)nb_allocate_zero_mem(size * sizeof(double));
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
		nb_free_mem(name);
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
		char* name = nb_allocate_mem(info[4]*sizeof(char));

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
				_A = nb_allocate_zero_mem(POW2(N) * sizeof(double));
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
		nb_free_mem(name);
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
	double* At = nb_allocate_zero_mem(POW2(N) * sizeof(double));
	uint32_t i, j;
	for(i=0; i < N; i++)
		for(j=0; j < N; j++)
			At[i*N + j] = A[j*N + i];

	/* Write matrix values */
	fwrite(At, sizeof(double), N*N, pfile);
	nb_free_mem(At);
	/* Close file */
	fclose(pfile);
}
