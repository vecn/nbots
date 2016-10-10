#ifndef __NB_MATH_BOT_NONPOLYNOMIAL_INTERPOLATOR_H__
#define __NB_MATH_BOT_NONPOLYNOMIAL_INTERPOLATOR_H__

void nb_nonpolynomial_eval(uint32_t N, uint8_t dim, const double *ni,
			   /* NULL for all ri = 1*/ const double *ri,
			   const double *x, double *eval);

void nb_nonpolynomial_eval_grad(uint32_t N, uint8_t dim,
				const double *ni,
				/* NULL for all ri = 1*/
				const double *ri, const double *x,
				double *eval);

#endif
