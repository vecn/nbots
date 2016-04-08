/******************************************************************************
 *   Sparse Bot: Linear Algebra for sparse and symmetric matrices.            *
 *   2011-2015 Victor Eduardo Cardoso Nungaray                                *
 *   Twitter: @victore_cardoso                                                *
 *   email: victorc@cimat.mx                                                  *
 ******************************************************************************/

/**
 * @file sparse_bot.h
 * @brief The sparse bot is a set of solvers for symmetric and sparse matrices.
 * There are fast matrix computations.
 * @author Victor Eduardo Cardoso Nungaray
 * @n victorc@@cimat.mx
 * @n <a href="https://twitter.com/victore_cardoso"> @@victore_cardoso </a>
 * @date 10 August 2015
 */

#ifndef __NB_EIGEN_BOT_SPARSE_H__
#define __NB_EIGEN_BOT_SPARSE_H__

#include <stdbool.h>
#include <stdint.h>
#include "nb/graph_bot.h"

  /**
   * @brief LU decomposition.
   */
#define NB_SOLVER_LUD (1) 

  /**
   * @brief Cholesky decomposition.
   */
#define NB_SOLVER_CHK (2)

  /**
   * @brief Conjugate Gradient preconditioned with Jacobi
   */
#define NB_SOLVER_CGJ (3)


/**
 * @brief Sparse matrix in "Compress Row Storage" format.
 */
typedef struct vcn_sparse_s vcn_sparse_t;
  
/**
 * @brief Create a sparse matrix from a graph.
 * @param[in] graph Create sparse matrix from a graph.
 * @param[in] perm Permutation (labeling) of the graph.
 * @param[in] vars_per_node Variables of the matrix per node in the graph.
 * @return Sparse matrix allocated or NULL if something goes wrong.
 */
vcn_sparse_t* vcn_sparse_create(const vcn_graph_t *const graph,
				const uint32_t *const perm,
				uint32_t vars_per_node);

vcn_sparse_t* vcn_sparse_clone(vcn_sparse_t* A);
void vcn_sparse_save(const vcn_sparse_t *const A, const char* filename);
void vcn_sparse_destroy(vcn_sparse_t* A);
void vcn_sparse_reset(vcn_sparse_t *A);
void vcn_sparse_set(vcn_sparse_t *A, uint32_t i, uint32_t j, double value);
void vcn_sparse_set_identity_row(vcn_sparse_t* A, uint32_t row);
void vcn_sparse_make_diagonal(vcn_sparse_t* A, double diag_val);
double vcn_sparse_get(const vcn_sparse_t *const A, uint32_t i, uint32_t j);
double vcn_sparse_get_and_set(vcn_sparse_t *A, uint32_t i, uint32_t j, double value);
bool vcn_sparse_is_non_zero(const vcn_sparse_t *const A, uint32_t i, uint32_t j);
uint32_t vcn_sparse_memory_used(const vcn_sparse_t *const A);
void vcn_sparse_add(vcn_sparse_t *A, uint32_t i, uint32_t j, double value);
void vcn_sparse_scale(vcn_sparse_t *A, double factor);
void vcn_sparse_transpose(vcn_sparse_t *A, vcn_sparse_t *_At);
uint32_t vcn_sparse_get_size(const vcn_sparse_t *const A);
uint32_t vcn_sparse_get_nnz(const vcn_sparse_t *const A);

void vcn_sparse_multiply_scalar(vcn_sparse_t* A, double scalar,
				uint32_t omp_parallel_threads);
void vcn_sparse_multiply_vector(vcn_sparse_t* A, double* in, double* out,
				uint32_t omp_parallel_threads);
  
/* Functions to explore sparse matrix */
int vcn_sparse_spy_plot_as_png(const vcn_sparse_t *const A,
			       const char* url, uint32_t size,
			       bool enable_zeros_allocated,
			       bool enable_color);

/* Set Dirichlet condition in the system */
void vcn_sparse_set_Dirichlet_condition(vcn_sparse_t* A, double* RHS,
					uint32_t idx, double value);

/* Iterative solvers for sparse matrix */
int vcn_sparse_solve_Gauss_Seidel
(const vcn_sparse_t *const A, 
 const double *const b,
 double *_x,                /* Out */
 uint32_t max_iter, double tolerance,
 uint32_t* niter_performed,     /* Out (NULL if not required) */
 double* tolerance_reached, /* Out (NULL if not required) */
 uint32_t omp_parallel_threads);

int vcn_sparse_solve_conjugate_gradient
(const vcn_sparse_t *const A, 
 const double *const b, 
 double *_x,                /* Out */
 uint32_t max_iter, double tolerance,
 uint32_t* niter_performed,     /* Out (NULL if not required) */
 double* tolerance_reached, /* Out (NULL if not required) */
 uint32_t omp_parallel_threads);

int vcn_sparse_solve_CG_precond_Jacobi
(const vcn_sparse_t *const A, 
 const double *const b,
 double *_x,                /* Out */
 uint32_t max_iter, double tolerance,
 uint32_t* niter_performed,     /* Out (NULL if not required) */
 double* tolerance_reached, /* Out (NULL if not required) */
 uint32_t omp_parallel_threads);

int vcn_sparse_solve_CG_precond_Cholesky
(const vcn_sparse_t *const A,
 const double *const b, 
 double *_x,                /* Out */
 uint32_t ktrunc,
 uint32_t max_iter, double tolerance,
 uint32_t* niter_performed,     /* Out (NULL if not required) */
 double* tolerance_reached, /* Out (NULL if not required) */
 uint32_t omp_parallel_threads);

int vcn_sparse_solve_CG_precond_fsai
(const vcn_sparse_t *const A,
 const double *const b,
 double *_x,                /* Out */
 double threshold,
 uint32_t max_iter,	double tolerance,
 uint32_t* niter_performed,     /* Out (NULL if not required) */
 double* tolerance_reached, /* Out (NULL if not required) */
 uint32_t omp_parallel_threads);

/* Factorization solvers for sparse matrix */
vcn_sparse_t* vcn_sparse_create_permutation
(const vcn_sparse_t *const A,
 const uint32_t *const perm,
 const uint32_t *const iperm);

void vcn_sparse_fill_permutation(const vcn_sparse_t *const A, 
				 vcn_sparse_t* Ar,
				 const uint32_t *const perm,
				 const uint32_t *const iperm);

double* vcn_sparse_create_vector_permutation
(const double *const b,
 const uint32_t *const perm, 
 uint32_t N);

/**
 * @brief Allocate LU decomposition using symbolic Cholesky
 * @param[in] A Matrix to be decomposed.
 * @param[in] L Pointer to the Lower triangular matrix that will be allocated.
 * @param[in] U Pointer to the Upper triangular matrix that will be allocated.
 * @return Zero if success. The matrix will be initialized in zeros.
 */
int vcn_sparse_alloc_LU
(const vcn_sparse_t *const A, vcn_sparse_t** L, vcn_sparse_t** U);

/**
 * @brief LU decomposition (Doolitle).
 * L and U must be allocated previously using vcn_sparse_alloc_LU().
 * @param[out] L Lower triangular matrix resulting from the decomposition.
 * @param[out] U Upper triangular matrix resulting from the decomposition.
 * @param[in] Number of threads to be used in parallel task.
 */
void vcn_sparse_decompose_LU(const vcn_sparse_t *const Ar,
			     vcn_sparse_t *L, vcn_sparse_t* U,
			     uint32_t omp_parallel_threads);

/**
 * @brief Cholesky decomposition.
 * L and Lt must be allocated previously using vcn_sparse_alloc_LU().
 * @param[out] L Lower triangular matrix resulting from the decomposition.
 * @param[out] Lt Upper triangular matrix resulting from the decomposition.
 * @param[in] Number of threads to be used in parallel task.
 */
int vcn_sparse_decompose_Cholesky(const vcn_sparse_t *const Ar, 
				  vcn_sparse_t * L, vcn_sparse_t* Lt,
				  uint32_t omp_parallel_threads);

void vcn_sparse_solve_LU(const vcn_sparse_t *const L, 
			 const vcn_sparse_t *const U,
			 const double *const b,
			 double* _x);  /* Out */
  
/* Default values solvers */
int vcn_sparse_solve_Cholesky(const vcn_sparse_t *const A,
			      const double *const b,
			      double* x,  /* Out */
			      uint32_t omp_parallel_threads);

int vcn_sparse_solve_using_LU(const vcn_sparse_t *const A,
			      const double *const b,
			      double* x,  /* Out */
			      uint32_t omp_parallel_threads);

/* Triangular solvers for sparse matrix */
void vcn_sparse_forward_solve(const vcn_sparse_t *const L,
			      const double *const b, 
			      double* _x /* Out */);
void vcn_sparse_backward_solve(const vcn_sparse_t *const U,
			       const double *const b,
			       double* _x /* Out */);

/* Eigenvalues and eigendoubles */
void vcn_sparse_eigen_power(const vcn_sparse_t *const A, int h,
			    double **_eigenvecs,/* Out */
			    double* _eigenvals, /* Out */
			    int* it,            /* Out */
			    double tolerance,
			    uint32_t omp_parallel_threads);
int vcn_sparse_eigen_ipower(const vcn_sparse_t *const A,
			    int solver_type,
			    int h, double mu,
			    double **_eigenvecs,/* Out */
			    double* _eigenvals, /* Out */
			    int* it,            /* Out */
			    double tolerance,
			    uint32_t omp_parallel_threads);
void vcn_sparse_eigen_lanczos(const vcn_sparse_t *const A,
			      double *_eigenmax,/* Out */ 
			      double* _eigenmin,/* Out */
			      int* it,          /* Out */
			      double tolerance,
			      uint32_t omp_parallel_threads);
void vcn_sparse_eigen_givens(const double* const main_diag, 
			     const double* const uplw_diag,
			     int i, double *_eigenvalue,
			     double tolerance,
			     uint32_t N);

/* Full matrix solvers */
void vcn_matrix_2X2_inverse(const double *const A,
			    double* A_inv);
  
void vcn_matrix_2X2_eigen(const double *const A,
			  /* Return: A decomposed into P Lambda P' */
			  double* Lambda,  /* Output (diag) */
			  double* P,       /* Output */
			  double tolerance);
double vcn_matrix_2X2_det(double *A);

double vcn_matrix_3X3_det(double *A);

void vcn_matrix_2X2_inverse_destructive(double *A);

void vcn_matrix_3X3_inverse_destructive(double *A);

int vcn_matrix_cholesky_decomposition
(const double *const A,
 double* _LplusLt,      /* Out */
 uint32_t N);
void vcn_matrix_cholesky_solve
(const double *const LplusLt,
 const double *const b, 
 double* _x,            /* Out */
 uint32_t N);
double vcn_matrix_cond1(const double *const A, int N);
double vcn_matrix_cond2(const double *const A, int N);
void vcn_matrix_qr_decomposition(double * A, /* Overwritten */
				 int N, 
				 double *c, double *d, int *sing);
void vcn_matrix_qr_solve(const double *const A,
			 int N, double *c, double *d,
			 double *b /* Solution overwritten */);
void vcn_matrix_svd_decomposition(double *A, /* Overwritten with U */
				  double *w, double *V, 
				  int N, int M);
void vcn_matrix_svd_solve(const double *const U,
			  const double *const w, 
			  const double *const V, 
			  double *x,             /* Out */
			  const double *const b, int N, int M);

/* Triangular solvers for complete matrix */
void vcn_matrix_forward_solve(const double *const L,
			      const double *const b,
			      double* _x, uint32_t N);
void vcn_matrix_backward_solve(const double *const U,
			       const double *const b,
			       double* _x, uint32_t N);

/* Functions to load and write Octave files (MAT-File 4 format) */
void vcn_sparse_read_mat4(vcn_sparse_t *A,
			  const char *url, char *label);
void vcn_sparse_save_mat4(const vcn_sparse_t *const A,
			  const char *url, char *label);

void vcn_mat4_printf(const char* url);
short vcn_mat4_exist(const char* url, char* label);
void vcn_mat4_clear(const char* url);
void vcn_mat4_read_vec(const char *url, char *label, double* _x);
void vcn_mat4_save_vec(const char *url, char *label,
		       const double *const  x, uint32_t N);
void vcn_mat4_read_mtx(const char *url, char *label, double* _A);
void vcn_mat4_save_mtx(const char *url, char *label,
		       const double *const A, uint32_t N);

#endif
