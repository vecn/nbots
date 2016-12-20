#ifndef __NB_SOLVER_BOT_VECTOR_H__
#define __NB_SOLVER_BOT_VECTOR_H__

#include <stdint.h>

double nb_vector_get_norm(const double *vec, uint32_t N);

void nb_vector_permutation(uint32_t N, const double *v,
			   const uint32_t *perm, double *vp);
#endif
