#ifndef __NB_PDE_BOT_FINITE_DIFFERENCES_GENERIC_ODE_H__
#define __NB_PDE_BOT_FINITE_DIFFERENCES_GENERIC_ODE_H__

#include <stdint.h>

/* -k f'' + a f' + b f = c */
void nb_fd_solve_ode_h2(uint32_t N, double L,
			double (*k)(double x), /* Diffusion (can be NULL) */
			double (*a)(double x), /* Convection (can be NULL) */
			double (*b)(double x), /* Absorption (can be NULL) */
			double (*c)(double x), /* Source (can be NULL) */
			double f0, double fL,  /* Boundary conditions */
			double *f_out);

#endif
