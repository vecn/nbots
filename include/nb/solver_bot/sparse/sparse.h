#ifndef __NB_SOLVER_BOT_SPARSE_SPARSE_H__
#define __NB_SOLVER_BOT_SPARSE_SPARSE_H__

#include <stdbool.h>
#include <stdint.h>

#include "nb/graph_bot.h"

typedef enum {
	NB_SOLVER_LUD,
	NB_SOLVER_CHK,
	NB_SOLVER_CGJ
} nb_solver_t;


/**
 * @brief Sparse matrix in "Compress Row Storage" format.
 */
typedef struct nb_sparse_s nb_sparse_t;
  
/**
 * @brief Create a sparse matrix from a graph.
 * @param[in] graph Create sparse matrix from a graph.
 * @param[in] perm Permutation (labeling) of the graph.
 * @param[in] vars_per_node Variables of the matrix per node in the graph.
 * @return Sparse matrix allocated or NULL if something goes wrong.
 */
nb_sparse_t* nb_sparse_create(const nb_graph_t *const graph,
			      const uint32_t *const perm,
			      uint32_t vars_per_node);

nb_sparse_t* nb_sparse_clone(nb_sparse_t* A);
void nb_sparse_save(const nb_sparse_t *const A, const char* filename);
void nb_sparse_destroy(nb_sparse_t* A);
void nb_sparse_reset(nb_sparse_t *A);
void nb_sparse_set(nb_sparse_t *A, uint32_t i, uint32_t j, double value);
void nb_sparse_set_identity_row(nb_sparse_t* A, uint32_t row);
void nb_sparse_make_diagonal(nb_sparse_t* A, double diag_val);
double nb_sparse_get(const nb_sparse_t *const A, uint32_t i, uint32_t j);
double nb_sparse_get_and_set(nb_sparse_t *A, uint32_t i, uint32_t j, double value);
bool nb_sparse_is_non_zero(const nb_sparse_t *const A, uint32_t i, uint32_t j);
uint32_t nb_sparse_memory_used(const nb_sparse_t *const A);
void nb_sparse_add(nb_sparse_t *A, uint32_t i, uint32_t j, double value);
void nb_sparse_scale(nb_sparse_t *A, double factor);
void nb_sparse_get_transpose(const nb_sparse_t *A, nb_sparse_t *_At);
void nb_sparse_transpose(nb_sparse_t *A);
uint32_t nb_sparse_get_size(const nb_sparse_t *const A);
uint32_t nb_sparse_get_nnz(const nb_sparse_t *const A);
double nb_sparse_get_asym(const nb_sparse_t *const A);
double nb_sparse_get_frobenius_norm(const nb_sparse_t *const A);

void nb_sparse_multiply_scalar(nb_sparse_t* A, double scalar,
				uint32_t omp_parallel_threads);
void nb_sparse_multiply_vector(const nb_sparse_t* A, const double* in,
				double* out, uint32_t omp_parallel_threads);
  
void nb_sparse_get_graph(const nb_sparse_t* A, nb_graph_t *graph);

/* Set Dirichlet condition in the system */
void nb_sparse_set_Dirichlet_condition(nb_sparse_t* A, double* RHS,
					uint32_t idx, double value);

/* Iterative solvers for sparse matrix */
int nb_sparse_solve_Gauss_Seidel
(const nb_sparse_t *const A, 
 const double *const b,
 double *_x,                /* Out */
 uint32_t max_iter, double tolerance,
 uint32_t* niter_performed,     /* Out (NULL if not required) */
 double* tolerance_reached, /* Out (NULL if not required) */
 uint32_t omp_parallel_threads);

int nb_sparse_solve_conjugate_gradient
(const nb_sparse_t *const A, 
 const double *const b, 
 double *_x,                /* Out */
 uint32_t max_iter, double tolerance,
 uint32_t* niter_performed,     /* Out (NULL if not required) */
 double* tolerance_reached, /* Out (NULL if not required) */
 uint32_t omp_parallel_threads);

int nb_sparse_solve_CG_precond_Jacobi
(const nb_sparse_t *const A, 
 const double *const b,
 double *_x,                /* Out */
 uint32_t max_iter, double tolerance,
 uint32_t* niter_performed,     /* Out (NULL if not required) */
 double* tolerance_reached, /* Out (NULL if not required) */
 uint32_t omp_parallel_threads);

int nb_sparse_solve_CG_precond_Cholesky
(const nb_sparse_t *const A,
 const double *const b, 
 double *_x,                /* Out */
 uint32_t ktrunc,
 uint32_t max_iter, double tolerance,
 uint32_t* niter_performed,     /* Out (NULL if not required) */
 double* tolerance_reached, /* Out (NULL if not required) */
 uint32_t omp_parallel_threads);

int nb_sparse_solve_CG_precond_fsai
(const nb_sparse_t *const A,
 const double *const b,
 double *_x,                /* Out */
 double threshold,
 uint32_t max_iter,	double tolerance,
 uint32_t* niter_performed,     /* Out (NULL if not required) */
 double* tolerance_reached, /* Out (NULL if not required) */
 uint32_t omp_parallel_threads);

/* Factorization solvers for sparse matrix */
nb_sparse_t* nb_sparse_create_permutation
(const nb_sparse_t *const A,
 const uint32_t *const perm,
 const uint32_t *const iperm);

void nb_sparse_fill_permutation(const nb_sparse_t *const A, 
				 nb_sparse_t* Ar,
				 const uint32_t *const perm,
				 const uint32_t *const iperm);

double* nb_sparse_create_vector_permutation
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
int nb_sparse_alloc_LU
(const nb_sparse_t *const A, nb_sparse_t** L, nb_sparse_t** U);

/**
 * @brief LU decomposition (Doolitle).
 * L and U must be allocated previously using nb_sparse_alloc_LU().
 * @param[out] L Lower triangular matrix resulting from the decomposition.
 * @param[out] U Upper triangular matrix resulting from the decomposition.
 * @param[in] Number of threads to be used in parallel task.
 */
void nb_sparse_decompose_LU(const nb_sparse_t *const Ar,
			     nb_sparse_t *L, nb_sparse_t* U,
			     uint32_t omp_parallel_threads);

/**
 * @brief Cholesky decomposition.
 * L and Lt must be allocated previously using nb_sparse_alloc_LU().
 * @param[out] L Lower triangular matrix resulting from the decomposition.
 * @param[out] Lt Upper triangular matrix resulting from the decomposition.
 * @param[in] Number of threads to be used in parallel task.
 */
int nb_sparse_decompose_Cholesky(const nb_sparse_t *const Ar, 
				  nb_sparse_t * L, nb_sparse_t* Lt,
				  uint32_t omp_parallel_threads);

void nb_sparse_solve_LU(const nb_sparse_t *const L, 
			 const nb_sparse_t *const U,
			 const double *const b,
			 double* _x);  /* Out */
  
/* Default values solvers */
int nb_sparse_solve_Cholesky(const nb_sparse_t *const A,
			      const double *const b,
			      double* x,  /* Out */
			      uint32_t omp_parallel_threads);

int nb_sparse_solve_using_LU(const nb_sparse_t *const A,
			      const double *const b,
			      double* x,  /* Out */
			      uint32_t omp_parallel_threads);

/* Triangular solvers for sparse matrix */
void nb_sparse_forward_solve(const nb_sparse_t *const L,
			      const double *const b, 
			      double* _x /* Out */);
void nb_sparse_backward_solve(const nb_sparse_t *const U,
			       const double *const b,
			       double* _x /* Out */);

/* Eigenvalues and eigendoubles */
void nb_sparse_eigen_power(const nb_sparse_t *const A, int h,
			    double **_eigenvecs,/* Out */
			    double* _eigenvals, /* Out */
			    int* it,            /* Out */
			    double tolerance,
			    uint32_t omp_parallel_threads);
int nb_sparse_eigen_ipower(const nb_sparse_t *const A,
			    nb_solver_t solver,
			    int h, double mu,
			    double **_eigenvecs,/* Out */
			    double* _eigenvals, /* Out */
			    int* it,            /* Out */
			    double tolerance,
			    uint32_t omp_parallel_threads);
void nb_sparse_eigen_lanczos(const nb_sparse_t *const A,
			      double *_eigenmax,/* Out */ 
			      double* _eigenmin,/* Out */
			      int* it,          /* Out */
			      double tolerance,
			      uint32_t omp_parallel_threads);
void nb_sparse_eigen_givens(const double* const main_diag, 
			     const double* const uplw_diag,
			     int i, double *_eigenvalue,
			     double tolerance,
			     uint32_t N);
#endif
