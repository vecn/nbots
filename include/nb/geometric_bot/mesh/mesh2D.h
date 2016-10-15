#ifndef __NB_GEOMETRIC_BOT_MESH_MESH2D_H__
#define __NB_GEOMETRIC_BOT_MESH_MESH2D_H__

#include <stdbool.h>
#include <stdint.h>

#include "nb/geometric_bot/model/model2D.h"
#include "nb/graph_bot.h"

/**
 * @brief Set as min_angle constrain in order to generate
 * a refined mesh with the maximum quality predicted by theory.
 * Minimum angle bound equivalent to 26.45 degrees.
 */
#define NB_MESH_MAX_ANGLE (0.46163958715250017309)

enum {
	NB_MESH_SIZE_CONSTRAINT_MAX_VTX,
	NB_MESH_SIZE_CONSTRAINT_MAX_TRG
};

enum {
	NB_MESH_GEOM_CONSTRAINT_MIN_ANGLE,
	NB_MESH_GEOM_CONSTRAINT_MAX_EDGE_LENGTH,
	NB_MESH_GEOM_CONSTRAINT_MAX_SUBSGM_LENGTH
};

enum {
	NB_MESH_TASK_AFTER_INSERT_TRG = 0,
	NB_MESH_TASK_AFTER_INSERT_VTX
};

enum {
	NB_MESH_REFINE_RUPPERT,
	NB_MESH_REFINE_CHEW,
	NB_MESH_REFINE_DEFAULT
};

/**
 * @brief Write-only mesh structure used to create and modify meshes.
 * This mesh is based on a Delaunay triangulation.
 */
typedef struct nb_tessellator2D__s nb_tessellator2D__t;

uint32_t nb_tessellator2D__get_memsize(void);
void nb_tessellator2D__init(nb_tessellator2D__t *mesh);
void nb_tessellator2D__finish(nb_tessellator2D__t *mesh);

nb_tessellator2D__t* nb_tessellator2D__create(void);
void nb_tessellator2D__clear(nb_tessellator2D__t* mesh);
void nb_tessellator2D__destroy(nb_tessellator2D__t* mesh);
void nb_tessellator2D__set_task(nb_tessellator2D__t *mesh, int type,
		       void (*task)(const nb_tessellator2D__t *const));

/**
 * @param[in] max_vtx Maximum number of vertices, which should be 
 * superior to the number of input vertices. Set zero for an undefined
 * maximum number of vertices (keeping guarantees).
 * The use of a max number of vertices could truncates the Delaunay 
 * refinment process if it runs-out of vertices. 
 * Use this feature only if it is strictly necessary.
 *
 * @param[in] max_trg Maximum number of triangles. Set zero for an undefined
 * maximum number of triangles (keeping guarantees).
 * The use of a max number of triangles could truncates the Delaunay 
 * refinment process if it runs-out of triangles.
 * Use this feature only if it is strictly necessary.
 */
void nb_tessellator2D__set_size_constraint(nb_tessellator2D__t *mesh, int type,
				  uint32_t value);
void nb_tessellator2D__unset_size_constraint(nb_tessellator2D__t *mesh, int type);
uint32_t nb_tessellator2D__get_size_constraint(const nb_tessellator2D__t *mesh, int type);

/**
 * @param[in] min_angle Minimum angle allowed in the triangulation (in 
 * radians).
 * This angle must be in the range of 0 and <b>NB_MESH_MAX_ANGLE</b>,
 * which corresponds to 26.45 degrees (0.4616 radians approx).
 */
void nb_tessellator2D__set_geometric_constraint(nb_tessellator2D__t *mesh, int type,
				       double value);
void nb_tessellator2D__unset_geometric_constraint(nb_tessellator2D__t *mesh, int type);
double nb_tessellator2D__get_geometric_constraint(const nb_tessellator2D__t *mesh,
					int type);
/**
 * @param[in] density Function to control the density of the triangulation,
 * which must be greater than zero for all points.
 * @n
 * Small <b>local features</b> could produce finer meshes than those 
 * indicated by the density function (because the minimum angle constraintt). 
 * @n
 * The expected size of an edge is the inverse of the density.
 * @n
 * The following values trigger a default behaviour:
 * + <b>NULL (or 0)</b>: No density defined, build the mesh with the minimum 
 *   number of triangles for a given minimum angle.
 *   <b>density_data</b> must be <b>NULL</b>.
 * + <b>NB_DENSITY_CDT</b>: Used to build the Constrainted Delaunay
 *   Triangulation (CDT) inside the domain. The CDT is the triangulation which
 *   maximizes the minimum angle of all the triangles, using only the input
 *   vertices.
 * + <b>NB_DENSITY_MAX</b>: Used to define maximum size of edges and input
 *   subsegments.
 *   <b>density_data</b> must be an array <b>'double max[2]'</b>, which
 *   defines the maximum edge length into max[0] and the maximum subsegment
 *   length into max[1]. To have any effect, the value max[1] must be smaller
 *   than max[0] because an input segment is also an edge (a greater value has
 *   not effect).
 *   If some max value equals zero then it is not considered, hence you can
 *   constraint only the subsegments size.
 * + <b>NB_DENSITY_IMG</b>: The density is given by an image.
 *   <b>density_data</b> must be a pointer to <b>nb_density_img_t</b>.
 *
 * @param[in] density_data Data used to estimate the density. It is given
 * to the density function without alterations.
 * Set NULL if not required.
 */
void nb_tessellator2D__set_density(nb_tessellator2D__t* mesh,
			  double (*density)(const double x[2],
					    const void *data),
			  const void *density_data);
void nb_tessellator2D__unset_density(nb_tessellator2D__t* mesh);
void nb_tessellator2D__set_refiner(nb_tessellator2D__t *mesh, int type);
int nb_tessellator2D__get_refiner(const nb_tessellator2D__t *const mesh);
bool nb_tessellator2D__is_empty(const nb_tessellator2D__t *const mesh);

/**
 * @brief Create an identical copy of the mesh.
 * @param[in] mesh Mesh to be cloned.
 * @return Cloned mesh if success, NULL if something goes wrong.
 */
nb_tessellator2D__t* nb_tessellator2D__clone(const nb_tessellator2D__t* const mesh);
  
/**
 * @brief Creates a mesh from a Planar Straight Line Graph (PSLG, encoded
 * in the model) as exposed by Schewchuck:
 *    + Calculates the Constrainted Delauany Triangulation from the PSLG.
 *    + Remove triangles from holes and concavities.
 *    + Delaunay Refinement (Circumcenter insertion of poor-quaility
 *      triangles). The internal routine is mesh_refine(), which
 *      is size optimal for a given minimum angle (please read the description
 *      of this routine).
 *
 * @see J. R. Shewchuk.
 * <b> Reprint of: Delaunay refinement algorithms for triangular mesh
 * generation.</b> Computational Geometry 47 (2014), pages 741-778.
 * @see G. L. Miller, S. E. Pav and N. J. Walkington.
 * <b> When and why Ruppert's algorithm works.</b> 2003
 *
 * @param[in, out] mesh Structure to store the mesh.
 * @param[in] model Structure which contains the PSLG to be triangulated;
 * the model also contains the location of the points inside the holes,
 * at least one point for each hole.
 */
void nb_tessellator2D__generate_from_model(nb_tessellator2D__t *mesh,
				  const nb_model_t *const model);

void nb_tessellator2D__get_simplest_from_model(nb_tessellator2D__t *mesh,
				      const nb_model_t *const  model);

/**
 * @brief Check if the vertex lies inside the mesh.
 * @param[in] mesh Discretization of the domain.
 * @param[in] vtx Vertex which is checked to be inside the domain.
 * @return <b>true</b> if the vertes lies inside the mesh 
 * (<b>false</b> if the vertex lies outside the mesh).
 */
bool nb_tessellator2D__is_vtx_inside(const nb_tessellator2D__t *const mesh,
			    const double *const vtx);

/**
 * @brief Produces a refined Delaunay triangulation. The algorithm have
 * the folowing properties:
 * + The minimal angle in the output is greater or equal to 26.45
 * degrees if the minimum input angle is greater than 36.53 degrees,
 * otherwise it produces the least quantity of poor-quality triangles.
 * + Termination is proven.
 * + Size optimality is proven.
 * + Good grading guarantees.
 *
 * @see J. Ruppert.
 * <b> A Delaunay Refinement Algorithm for Quality 2-Dimensional
 * Mesh Generation.</b> Journal of Algorithms 18 (1995), pages 548-585.
 * @see J. R. Shewchuk.
 * <b> Reprint of: Delaunay refinement algorithms for triangular mesh
 * generation.</b> Computational Geometry 47 (2014), pages 741-778.
 * @see G. L. Miller, S. E. Pav and N. J. Walkington.
 * <b> When and why Ruppert's algorithm works.</b> 2003
 *
 * @param[in] mesh The triangulation which is going to be refined, it
 * is assumed to be a Constrainted Delaunay Triangulation.
 */
void nb_tessellator2D__refine(nb_tessellator2D__t *mesh);
bool nb_tessellator2D__insert_vtx(nb_tessellator2D__t *mesh, const double vertex[2]);
void nb_tessellator2D__get_vertices(nb_tessellator2D__t* mesh, double* vertices);
uint32_t nb_tessellator2D__get_N_vtx(const nb_tessellator2D__t *const mesh);
uint32_t nb_tessellator2D__get_N_trg(const nb_tessellator2D__t *const mesh);
uint32_t nb_tessellator2D__get_N_edg(const nb_tessellator2D__t *const mesh);
double nb_tessellator2D__get_area(const nb_tessellator2D__t *const mesh);

#endif
