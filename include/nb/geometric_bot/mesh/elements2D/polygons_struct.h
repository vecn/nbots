#ifndef __VCN_MESH2D_POLYGONS_STRUCT_H__
#define __VCN_MESH2D_POLYGONS_STRUCT_H__

/**
 * @brief Read-only mesh structure storing the Voronoi graph, which
 * is a dual tesselation of the Delaunay triangulation. The Voronoi graph
 * is a mesh conformed by convex polygons, with the property that the closest
 * centroid of any point in the graph is the centroid of the polygon where
 * such a point is located.
 */
typedef struct vcn_mshpoly_s {
	uint32_t N_vertices;
	double* vertices;   /* Vertices coordinates concatenated */
	uint32_t N_polygons;
	double* centroids;  /* Centroids coordinates concatenated */
	uint32_t* N_vertices_forming_polygons;
	uint32_t** vertices_forming_polygons; /* Connectivity matrix */
	uint32_t* N_adjacencies;
	uint32_t** adjacencies;
} vcn_mshpoly_t;


#endif
