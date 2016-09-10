#ifndef __NB_PDE_BOT_FRECHET_DERIVATIVE_H__
#define __NB_PDE_BOT_FRECHET_DERIVATIVE_H__

void nb_pde_get_frechet_derivative(uint8_t N, uint8_t dim_x, uint8_t dim_f,
				   const double *x, const double *f,
				   const double *ni, const double *fi,
				   double *Df);

#endif
