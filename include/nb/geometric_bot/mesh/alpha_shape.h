#ifndef __NB_GEOMETRIC_BOT_MESH_ALPHA_SHAPE_H__
#define __NB_GEOMETRIC_BOT_MESH_ALPHA_SHAPE_H__

#include "nb/geometric_bot/mesh/tessellator2D.h"

void nb_tessellator2D_get_alpha_complex(nb_tessellator2D_t *mesh,
					uint32_t N_vertices,
					const double *const vertices,
					double alpha);
/**
 * @brief Smallest non-singular
 */
double nb_tessellator2D_get_smallest_ns_alpha_complex(nb_tessellator2D_t *mesh,
					     uint32_t N_vertices,
					     const double *vertices,
					     double alpha_factor);

#endif
