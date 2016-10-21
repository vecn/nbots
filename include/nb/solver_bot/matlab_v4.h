#ifndef __NB_SOLVER_BOT_MATLAB_V4_H__
#define __NB_SOLVER_BOT_MATLAB_V4_H__

#include "nb/solver_bot/sparse/sparse.h"

#include <stdint.h>

void nb_sparse_read_mat4(nb_sparse_t *A, const char *url,
			 const char *label);
void nb_sparse_save_mat4(const nb_sparse_t *const A,
			  const char *url, char *label);

void nb_mat4_printf(const char* url);
short nb_mat4_exist(const char* url, char* label);
void nb_mat4_clear(const char* url);
void nb_mat4_read_vec(const char *url, char *label, double* _x);
void nb_mat4_save_vec(const char *url, char *label,
		       const double *const  x, uint32_t N);
void nb_mat4_read_mtx(const char *url, char *label, double* _A);
void nb_mat4_save_mtx(const char *url, char *label,
		       const double *const A, uint32_t N);

#endif
