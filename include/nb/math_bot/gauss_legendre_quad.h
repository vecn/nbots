#ifndef __NB_MATH_BOT_GAUSS_LEGENDRE_QUAD_H__
#define __NB_MATH_BOT_GAUSS_LEGENDRE_QUAD_H__

#include <stdint.h>

typedef struct {
	uint8_t N;
	double* w;
	double* x;
} nb_glquadrature_t;

void nb_glquadrature_load(nb_glquadrature_t *glq, uint8_t N_points);

#endif
