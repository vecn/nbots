#ifndef __NB_PDE_BOT_ODE_SOLVERS_H__
#define __NB_PDE_BOT_ODE_SOLVERS_H__

#include <stdint.h>

/* -k(x) f'' + a(x) f' + b(x) f = c(x) */
void nb_fd_solve_ode_h2(uint32_t N, double L,
			double (*k)(double x), /* Diffusion (can be NULL) */
			double (*a)(double x), /* Convection (can be NULL) */
			double (*b)(double x), /* Absorption (can be NULL) */
			double (*c)(double x), /* Source (can be NULL) */
			double f0, double fL,  /* Boundary conditions */
			double *f_out);

/* -k(x) f'' + a(x) f' + b(x) f = c(x) */
void nb_fem_solve_ode_o1(uint32_t N, double L,
			 double (*k)(double x), /* Diffusion (can be NULL) */
			 double (*a)(double x), /* Convection (can be NULL) */
			 double (*b)(double x), /* Absorption (can be NULL) */
			 double (*c)(double x), /* Source (can be NULL) */
			 double f0, double fL,  /* Boundary conditions */
			 double *f_out);

/* -k(x) f'' + a(x) f' = c(x) */
void nb_analytic_fd_solve_diffusion_convection
			(uint32_t N, double L,
			 double (*k)(double x), /* Diffusion (can be NULL) */
			 double (*a)(double x), /* Convection (can be NULL) */
			 double (*c)(double x), /* Source (can be NULL) */
			 double f0, double fL,  /* Boundary conditions */
			 double *f_out);

/* -k(x) f'' + a(x) f' = c(x) */
void nb_analytic_fem_solve_diffusion_convection
			(uint32_t N, double L,
			 double (*k)(double x), /* Diffusion (can be NULL) */
			 double (*a)(double x), /* Convection (can be NULL) */
			 double (*c)(double x), /* Source (can be NULL) */
			 double f0, double fL,  /* Boundary conditions */
			 double *f_out);

#endif
