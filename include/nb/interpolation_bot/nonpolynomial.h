#ifndef __NB_INTERPOLATION_BOT_NONPOLYNOMIAL_H__
#define __NB_INTERPOLATION_BOT_NONPOLYNOMIAL_H__

void nb_nonpolynomial_simple_eval(uint8_t N, uint8_t dim, const double *ni,
				  const double *x, double *eval);
void nb_nonpolynomial_simple_eval_grad(uint8_t N, uint8_t dim,
				       const double *ni, const double *x,
				       double *eval);
void nb_nonpolynomial_thin_plate_eval(uint8_t N, uint8_t dim, const double *ni,
				      const double *x, double *eval);
void nb_nonpolynomial_thin_plate_eval_grad(uint8_t N, uint8_t dim,
					   const double *ni, const double *x,
					   double *eval);

void nb_nonpolynomial_custom_eval(uint8_t N, uint8_t dim, const double *ni,
				  /* NULL for all ri = 1*/ const double *ri,
				  const double *x, double *eval,
				  /* NULL for g(x) = x */ double (*g)(double));

void nb_nonpolynomial_custom_eval_grad(uint8_t N, uint8_t dim,
				       const double *ni,
				       /* NULL for all ri = 1*/
				       const double *ri,
				       const double *x, double *eval,
				       /* NULL for g(x) = x */
				       double (*g)(double),
				       /* NULL for g'(x) = 1 */
				       double (*dg)(double));

#endif
