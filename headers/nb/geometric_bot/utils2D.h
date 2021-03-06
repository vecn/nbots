#ifndef __NB_GEOMETRIC_BOT_UTILS_2D_H__
#define __NB_GEOMETRIC_BOT_UTILS_2D_H__

#include <stdbool.h>
#include <stdint.h>

#define NB_GEOMETRIC_TOL (1.5e-8)
#define NB_GEOMETRIC_TOL_POW2 (2.25e-16)

typedef enum {
	NB_INTERSECTED,
	NB_NOT_INTERSECTED,
	NB_PARALLEL,
	NB_INTERSECT_ON_A1,
	NB_INTERSECT_ON_A2,
	NB_INTERSECT_ON_B1,
	NB_INTERSECT_ON_B2
} nb_intersect_t;

void nb_utils2D_init(void);

double nb_utils2D_get_x_from_darray(const void *const vtx_ptr);
double nb_utils2D_get_y_from_darray(const void *const vtx_ptr);
void nb_utils2D_get_normal(const double x1[2],
			   const double x2[2],
			   double normal[2]);
double nb_utils2D_get_dist(const double p1[2], const double p2[2]);
double nb_utils2D_get_dist2(const double p1[2], const double p2[2]);
double nb_utils2D_get_delaunay_dist(const double p1[2],
				    const double p2[2],
				    const double p3[2]);
double nb_utils2D_get_2vec_angle(const double a[2],
				 const double b[2],
				 const double c[2]);
void nb_utils2D_get_enveloping_box
(uint32_t N_vertices,
 const void* vertices,
 uint16_t size_of_vtx_type,
 double (*get_x)(const void *const vtx),
 double (*get_y)(const void *const vtx),
 double box[4]);
void nb_utils2D_get_enveloping_box_from_subset
(uint32_t N_subset,
 const uint32_t subset[],
 const void* vertices,
 uint16_t size_of_vtx_type,
 double (*get_x)(const void *const vtx),
 double (*get_y)(const void *const vtx),
 double box[4]);
double nb_utils2D_orient(const double t1[2],
			 const double t2[2],
			 const double t3[2]);
bool nb_utils2D_is_in_half_side(const double v1[2],
				const double v2[2],
				const double v3[2]);
double nb_utils2D_get_trg_area(const double t1[2],
			       const double t2[2],
			       const double t3[2]);
double nb_utils2D_get_poly_area(const double *p[2], uint16_t N);
void nb_utils2D_trg_get_centroid(const double t1[2],
				 const double t2[2],
				 const double t3[2],
				 double centroid[2]);
double nb_utils2D_get_circumradius(const double t1[2],
				   const double t2[2],
				   const double t3[2]);
void nb_utils2D_get_circumcenter(const double t1[2],
				 const double t2[2],
				 const double t3[2],
				 double circumcenter[2]);
void nb_utils2D_get_circumcenter_from_sgm(const double s1[2],
					  const double s2[2],
					  double radii,
					  double circumcenter[2]);
double nb_utils2D_get_trg_min_angle(const double t1[2],
				    const double t2[2],
				    const double t3[2]);
double nb_utils2D_get_cr2se_ratio(const double t1[2],
				  const double t2[2],
				  const double t3[2]);
double nb_utils2D_get_trg_quality(const double t1[2],
				  const double t2[2],
				  const double t3[2]);
double nb_utils2D_get_min_trg_edge(const double t1[2],
				   const double t2[2],
				   const double t3[2]);
double nb_utils2D_get_max_trg_edge(const double t1[2],
				   const double t2[2],
				   const double t3[2]);
uint32_t nb_utils2D_get_closest_vtx(double x, double y,
				    uint32_t N_vertices,
				    const double vertices[]);
uint32_t nb_utils2D_get_closest_vtx_not_ignored
(double x, double y,
 uint32_t N_vertices,
 const double vertices[],
 uint32_t N_index_to_ignore,
 const uint32_t index_to_ignore[]);
void nb_utils2D_get_closest_pnt_to_sgm(const double s1[2],
				       const double s2[2],
				       const double p[2],
				       double closest_point[2]);
nb_intersect_t nb_utils2D_get_sgm_intersection(const double a1[2],
					       const double a2[2],
					       const double b1[2],
					       const double b2[2],
					       /* Output NULL if not required */
					       double intersection[2]);
bool nb_utils2D_are_sgm_intersected(const double a1[2], const double a2[2],
				    const double b1[2], const double b2[2],
				    /* Output NULL if not required */
				    double intersection[2]);
bool nb_utils2D_sgm_intersects_trg(const double t1[2],
				   const double t2[2],
				   const double t3[2],
				   const double s1[2],
				   const double s2[2]);

bool nb_utils2D_sgm_intersects_circle(const double circumcenter[2],
				      double radius,
				      const double s1[2],
				      const double s2[2]);

bool nb_utils2D_level_set_intersects_trg(double v1, double v2, double v3,
					 double level_set);
void nb_utils2D_get_trg_level_set_intersection(const double t1[2],
					       const double t2[2],
					       const double t3[2],
					       double v1, double v2, double v3,
					       double level_set,
					       double a[2], double b[2]);

bool nb_utils2D_pnt_lies_on_sgm(const double s1[2],
				const double s2[2],
				const double p[2]);

bool nb_utils2D_pnt_lies_in_trg(const double t1[2],
				const double t2[2],
				const double t3[2],
				const double p[2]);

bool nb_utils2D_pnt_lies_in_poly(int N, const double *poly,
				 const double p[2]);
bool nb_utils2D_pnt_lies_in_poly_bnd(int N, const double *poly,
				     const double p[2]);

bool nb_utils2D_pnt_lies_in_diametral_circle
(const double s1[2],
 const double s2[2],
 const double p[2]);

bool nb_utils2D_pnt_lies_strictly_in_diametral_circle
(const double s1[2],
 const double s2[2],
 const double p[2]);

bool nb_utils2D_pnt_lies_strictly_in_circumcircle(const double t1[2],
						  const double t2[2],
						  const double t3[2],
						  const double p[2]);

bool nb_utils2D_pnt_is_cocircular(const double t1[2],
				  const double t2[2],
				  const double t3[2],
				  const double p[2]);

bool nb_utils2D_pnt_lies_in_box(const double box[4],
				const double p[2]);

#endif
