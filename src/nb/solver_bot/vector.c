#include <stdint.h>
#include <math.h>

#include "nb/solver_bot/vector.h"

#define POW2(a) ((a)*(a))

double nb_vector_get_norm(const double* x, uint32_t N)
{
	double n = 0;
	for (uint32_t i = 0; i < N; i++)
		n += POW2(x[i]);
	return sqrt(n);
}
