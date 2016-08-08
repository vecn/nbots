#include <math.h>

#include "nb/pde_bot/solid_mechanics_formulas.h"

#define POW2(a) ((a)*(a))

double nb_pde_get_vm_stress(double sxx, double syy, double sxy)
{
	return sqrt(POW2(sxx) + POW2(syy) - sxx * syy + 3.0 * POW2(sxy));
}
