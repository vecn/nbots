#ifndef __NB_INTERPOLATION_BOT_NONPOLYNOMIAL_H__
#define __NB_INTERPOLATION_BOT_NONPOLYNOMIAL_H__

void nb_nonpolynomial_eval(uint8_t N, uint8_t dim, const double *ni,
			   /* NULL for all ri = 1*/ const double *ri,
			   const double *x, double c, double *eval);

void nb_nonpolynomial_eval_grad(uint8_t N, uint8_t dim,
				const double *ni,
				/* NULL for all ri = 1*/
				const double *ri, const double *x,
				double c, double *eval);

#endif
