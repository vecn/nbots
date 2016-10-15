#ifndef __NB_GEOMETRIC_BOT_MESH_ALPHA_SHAPE_H__
#define __NB_GEOMETRIC_BOT_MESH_ALPHA_SHAPE_H__

#include "nb/geometric_bot/mesh/mesh2D.h"

void nb_tessellator2D__get_alpha_complex(nb_tessellator2D__t *mesh, uint32_t N_vertices,
			       const double *const vertices,
			       double alpha);
/**
 * @brief Smallest non-singular
 */
double nb_tessellator2D__get_smallest_ns_alpha_complex(nb_tessellator2D__t *mesh,
					     uint32_t N_vertices,
					     const double *vertices,
					     double alpha_factor);

#endif
