#include <stdlib.h>

#include "nb/geometric_bot.h"

int nb_cvfa_interpolate_from_centroids_to_nodes
		(const nb_partition_t *const part,
		 uint32_t N_comp,
		 const double* cen_values,
		 double* nodal_values /* Output */)
{
	/* Use shape functions of Sukumar??? No support for 3D */
	/* Use rational polynomials?? */
	return 0;/* PENDING */
}
