#ifndef __VCN_GEOMETRIC_BOT_MESH_MESH2D_H__
#define __VCN_GEOMETRIC_BOT_MESH_MESH2D_H__

#include <stdbool.h>
#include <stdint.h>

#include "vcn/geometric_bot/model/model2D.h"
#include "vcn/graph_bot.h"

/**
 * @brief Set as min_angle constrain in order to generate
 * a refined mesh with the maximum quality predicted by theory.
 * Minimum angle bound equivalent to 26.45 degrees.
 */
#define VCN_MESH_MAX_ANGLE (0.46163958715250017309)

enum {
	VCN_MESH_SIZE_CONSTRAINT_MAX_VTX,
	VCN_MESH_SIZE_CONSTRAINT_MAX_TRG
};

enum {
	VCN_MESH_GEOM_CONSTRAINT_MIN_ANGLE,
	VCN_MESH_GEOM_CONSTRAINT_MAX_EDGE_LENGTH,
	VCN_MESH_GEOM_CONSTRAINT_MAX_SUBSGM_LENGTH
};

enum {
	VCN_MESH_TASK_AFTER_INSERT_TRG,
	VCN_MESH_TASK_AFTER_INSERT_VTX
};

enum {
	VCN_MESH_REFINE_RUPPERT,
	VCN_MESH_REFINE_CHEW,
	VCN_MESH_REFINE_DEFAULT
};

/**
 * @brief Write-only mesh structure used to create and modify meshes.
 * This mesh is based on a Delaunay triangulation.
 * In order to read the mesh, it must be "converted" into a read-only 
 * structure, such as vcn_msh3trg_t.
 */
typedef struct vcn_mesh_s vcn_mesh_t;

vcn_mesh_t* vcn_mesh_create(void);
void vcn_mesh_clear(vcn_mesh_t* mesh);
void vcn_mesh_destroy(vcn_mesh_t* mesh);
void vcn_mesh_set_task(vcn_mesh_t *mesh, int type,
		       void (*task)(const vcn_mesh_t *const));
void vcn_mesh_set_size_constraint(vcn_mesh_t *mesh, int type,
				 uint32_t value);
uint32_t vcn_mesh_get_size_constraint(const vcn_mesh_t *mesh, int type);
void vcn_mesh_set_geometric_constraint(vcn_mesh_t *mesh, int type,
				      double value);
double vcn_mesh_get_geometric_constraint(const vcn_mesh_t *mesh,
					int type);
void vcn_mesh_set_density(vcn_mesh_t* mesh,
			  double (*density)(const double x[2],
					    const void *data),
			  const void *density_data);
void vcn_mesh_set_refiner(vcn_mesh_t *mesh, int type);
int vcn_mesh_get_refiner(const vcn_mesh_t *const mesh);
bool vcn_mesh_is_empty(const vcn_mesh_t *const mesh);

/**
 * @brief Create an identical copy of the mesh.
 * @param[in] mesh Mesh to be cloned.
 * @return Cloned mesh if success, NULL if something goes wrong.
 */
vcn_mesh_t* vcn_mesh_clone(const vcn_mesh_t* const mesh);
  
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
void vcn_mesh_generate_from_model(vcn_mesh_t *mesh,
				  const vcn_model_t *const model);

/**
 * @brief Check if the vertex lies inside the mesh.
 * @param[in] mesh Discretization of the domain.
 * @param[in] vtx Vertex which is checked to be inside the domain.
 * @return <b>true</b> if the vertes lies inside the mesh 
 * (<b>false</b> if the vertex lies outside the mesh).
 */
bool vcn_mesh_is_vtx_inside(const vcn_mesh_t *const mesh,
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
 *
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
 *
 * @param[in] min_angle Minimum angle allowed in the triangulation (in 
 * radians).
 * This angle must be in the range of 0 and <b>VCN_MESH_MAX_ANGLE</b>,
 * which corresponds to 26.45 degrees (0.4616 radians approx).
 *
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
 * + <b>VCN_DENSITY_CDT</b>: Used to build the Constrainted Delaunay
 *   Triangulation (CDT) inside the domain. The CDT is the triangulation which
 *   maximizes the minimum angle of all the triangles, using only the input
 *   vertices.
 * + <b>VCN_DENSITY_MAX</b>: Used to define maximum size of edges and input
 *   subsegments.
 *   <b>density_data</b> must be an array <b>'double max[2]'</b>, which
 *   defines the maximum edge length into max[0] and the maximum subsegment
 *   length into max[1]. To have any effect, the value max[1] must be smaller
 *   than max[0] because an input segment is also an edge (a greater value has
 *   not effect).
 *   If some max value equals zero then it is not considered, hence you can
 *   constraint only the subsegments size.
 * + <b>VCN_DENSITY_IMG</b>: The density is given by an image.
 *   <b>density_data</b> must be a pointer to <b>vcn_density_img_t</b>.
 *
 * @param[in] density_data Data used to estimate the density. It is given
 * to the density function without alterations.
 * Set NULL if not required.
 */
void vcn_mesh_refine(vcn_mesh_t *mesh);
bool vcn_mesh_insert_vtx(vcn_mesh_t *mesh, const double vertex[2]);
void vcn_mesh_get_vertices(vcn_mesh_t* mesh, double* vertices);
uint32_t vcn_mesh_get_N_vtx(const vcn_mesh_t *const mesh);
uint32_t vcn_mesh_get_N_trg(const vcn_mesh_t *const mesh);
uint32_t vcn_mesh_get_N_edg(const vcn_mesh_t *const mesh);
double vcn_mesh_get_area(const vcn_mesh_t *const mesh);

#endif
