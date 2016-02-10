#ifndef __NB_GEOMETRIC_BOT_MESH_ELEMENTS2D_DISKS_H__
#define __NB_GEOMETRIC_BOT_MESH_ELEMENTS2D_DISKS_H__

#include <stdint.h>

#include "disks_struct.h"

#include "nb/geometric_bot/mesh/mesh2D.h"
#include "nb/graph_bot.h"

vcn_mshpack_t* vcn_mesh_get_mshpack
(const vcn_mesh_t *const mesh,
 bool include_adjacencies,
 uint32_t iterations,
 double overlapping_factor,  /* Overlapping percentage [0,1] */
 double porosity_factor,     /* Porosity percentage [0,1] */
 uint32_t* (*labeling)(const vcn_graph_t *const) /* NULL for an arbitrary labeling */);

void vcn_mshpack_clear_adj(vcn_mshpack_t* spack);

/*! \fn void vcn_mshpack_destroy(vcn_mshpack_t* spack)
 */
void vcn_mshpack_destroy(vcn_mshpack_t* spack);

#endif
