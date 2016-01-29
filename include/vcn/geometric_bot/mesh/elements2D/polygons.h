#ifndef __VCN_MESH2D_POLYGONS_H__
#define __VCN_MESH2D_POLYGONS_H__

#include <stdbool.h>
#include <stdint.h>

#include "polygons_struct.h"
#include "vcn/geometric_bot/mesh/mesh2D.h"
#include "vcn/graph_bot.h"

#ifdef __cplusplus
extern "C" {
#endif

  vcn_mshpoly_t* vcn_mesh_get_mshpoly
        (const vcn_mesh_t *const mesh,
	 bool include_adjacencies,
	 bool central_voronoi,
	 uint32_t central_voronoi_max_iter,
	 /* NULL for a constant density */
	 double (*central_voronoi_density)(double*),
	 uint32_t* (*labeling)(const vcn_graph_t *const));/* NULL for an arbitrary labeling */

  void vcn_mshpoly_destroy(vcn_mshpoly_t* voronoi);

#ifdef __cplusplus
}
#endif

#endif
