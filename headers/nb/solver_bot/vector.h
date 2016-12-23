#ifndef __NB_SOLVER_BOT_VECTOR_H__
#define __NB_SOLVER_BOT_VECTOR_H__

#include <stdint.h>

double nb_vector_get_norm(const double *vec, uint32_t N);

void nb_vector_permutation(uint32_t N, const double *v,
			   const uint32_t *perm, double *vp);
void nb_vector_sum(uint32_t N, double *a, const double *b);
void nb_vector_substract(uint32_t N, double *a, const double *b);
void nb_vector_substract_to(uint32_t N, double *a, const double *b);
void nb_vector_scale(uint32_t N, double *a, double factor);

#endif
