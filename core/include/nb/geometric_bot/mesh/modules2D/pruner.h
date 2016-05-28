#ifndef __NB_GEOMETRIC_BOT_MESH_MODULES2D_PRUNER_H__
#define __NB_GEOMETRIC_BOT_MESH_MODULES2D_PRUNER_H__

#include <stdint.h>
#include <stdbool.h>

/**
 * @brief Calculates the pseudo-centroid of the sub-areas in the mesh.
 * A sub-area is a portion of the domain which is bounded by input segments,
 * A single mesh could have multiple subareas.
 * The pseudo-centroid is the centroid of a triangle (inside the sub-area)
 * which is closest to the real centroid of the whole subarea.
 * The pseudo-centroid is guaranteed to be inside a concave polygon, such as
 * most of the sub-areas.
 * @param[in] mesh Mesh to analyse.
 * @param[out] N_centroids Stores the number of enveloped areas.
 * @return Return the array with the locations of the centroids concatenated.
 */  
double* vcn_mesh_get_centroids_of_subareas(const vcn_mesh_t *const mesh,
					   uint32_t* N_centroids);

/**
 * @brief A hole can be seen as an enveloped sub-area, which is a sub-area inside
 * another enveloping subarea. This routine calculates the location of the 
 * pseudo-centroids of enveloped areas, which could be used to plant holes.
 * @param[in] mesh Mesh to analyse.
 * @param[out] N_centroids Stores the number of enveloped areas.
 * @return Return the array with the locations of the centroids concatenated.
 */
double* vcn_mesh_get_centroids_of_enveloped_areas(const vcn_mesh_t *const mesh,
						  uint32_t* N_centroids);


/**
 * @brief Delete the triangles inside enveloped areas.
 * @param[in] mesh Mesh to be processed.
 * @param[out] area_removed NULL if not required.
 * Stores the total area of deleted triangles.
 * @return Area of the non-deleted mesh.
 */
double vcn_mesh_clear_enveloped_areas(vcn_mesh_t* mesh,
				      double* area_removed);

/**
 * @brief If the mesh is not a continuous domain, this function keep 
 * the biggest sub-area and delete the other ones.
 * @n<b>WARNING</b>: After this functions is used, some input segments and input
 * vertices will be disconnected.
 * @param[in] mesh Mesh to be processed.
 * @param[out] area_removed NULL if not required.
 * Stores the total area of deleted triangles.
 * @return Area of the non-deleted mesh.
 */
double vcn_mesh_keep_biggest_continuum_area(vcn_mesh_t* mesh,
					    double* area_removed);

/**
 * @brief Delete those input segments without adjacent triangles, which could arise
 * if the PSLG is not well defined or if one of the following functions has been used:
 * vcn_mesh_clear_enveloped_areas() or vcn_mesh_keep_biggest_isolated_area().
 * @param[in] mesh Mesh to be processed.
 * @return Number of subsegments deleted.
 */
uint32_t vcn_mesh_delete_isolated_segments(vcn_mesh_t *const mesh);

  
/**
 * @brief Delete those input segments inside the boundaries.
 * @param[in] mesh Mesh to be processed.
 * @return Number of subsegments deleted.
 */
uint32_t vcn_mesh_delete_internal_input_segments(vcn_mesh_t *const mesh);

/**
 * @brief Delete those input vertices disconnected from the mesh, which could arise
 * if the PSLG is not well defined or if the function
 * vcn_mesh_delete_isolated_segments() has been used.
 * @param[in] mesh Mesh to be processed.
 * @return Number of vertices deleted.
 */
uint32_t vcn_mesh_delete_isolated_vertices(vcn_mesh_t* mesh);

bool vcn_mesh_is_continuum(const vcn_mesh_t *mesh);

uint16_t vcn_mesh_get_N_subareas(const vcn_mesh_t *mesh);

uint16_t nb_mesh_get_subareas(const vcn_mesh_t *mesh, uint16_t *area_id);

uint16_t vcn_mesh_get_N_continuum_areas(const vcn_mesh_t *mesh);

#endif
