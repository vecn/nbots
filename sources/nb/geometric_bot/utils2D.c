#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "nb/math_bot.h"
#include "nb/container_bot.h"
#include "nb/geometric_bot/utils2D.h"

#include "tiny_libs/predicates.h"

#define INCIRCLE_TOLERANCE 1e-9

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define POW2(a) ((a)*(a))

static void init_robust_predicates_if_not(void);
static char init_id = 0;

static inline bool level_set_intersects_sgm(double v1, double v2,
					    double level_set);
static void get_sgm_level_set_intersection(const double s1[2],
					   const double s2[2],
					   double v1, double v2,
					   double level_set,
					   double p[2]);

static bool right_ray_is_intersected(const double v1[2], const double v2[2],
				     const double p[2]);

static double det_circumcircle(const double t1[2],
			       const double t2[2],
			       const double t3[2],
			       const double p[2]);

double nb_utils2D_get_x_from_darray(const void *const vtx_ptr)
{
	const double *const vtx = vtx_ptr;
	return vtx[0];
}

double nb_utils2D_get_y_from_darray(const void *const vtx_ptr)
{
	const double *const vtx = vtx_ptr;
	return vtx[1];
}

void nb_utils2D_get_normal(const double x1[2], const double x2[2],
			   double normal[2])
/* Return normal from x1 towards x2 */
{
	normal[0] = x2[0] - x1[0];
	normal[1] = x2[1] - x1[1];
	double normalizer = sqrt(POW2(normal[0]) +
				 POW2(normal[1]));
	normal[0] /= normalizer;
	normal[1] /= normalizer;
}

double nb_utils2D_get_dist(const double p1[2], const double p2[2])
{
	return nb_math_hypo(p2[0]-p1[0], p2[1]-p1[1]);
}

double nb_utils2D_get_dist2(const double p1[2], const double p2[2])
{
	return POW2(p2[0]-p1[0]) + POW2(p2[1]-p1[1]);
}

double nb_utils2D_get_delaunay_dist(const double p1[2], const double p2[2],
				     const double p3[2])
{
	/* The Delaunay distance is defined by:
	 *
	 *                         | r**2  if  cc in Halfspace(v1, v2, v3)
	 *       dd(v1, v2, v3) =  |-r**2  otherwise,
	 *
	 * where cc is the circumcenter and Halfspace(v1, v2, v3) is
	 * the side of the wall where the triangle lies.
	 */
	double cc[2];
	nb_utils2D_get_circumcenter(p1, p2, p3, cc);
	double circumradius = nb_utils2D_get_dist2(cc, p1);
	double sign = nb_utils2D_orient(p1, p2, cc);
	if(sign < 0.0)
		circumradius *= -1;
	return circumradius;
}

double nb_utils2D_get_2vec_angle(const double a[2],
				 const double b[2],
				 const double c[2])
{
	/* Angle between segments A and B, formed by vtx a, b and c
	 *           a       c
	 *            \     /
	 *           A \___/ B
	 *              \ /
	 *               b
	 */
	double A[2];
	A[0] = a[0] - b[0];
	A[1] = a[1] - b[1];

	double B[2];
	B[0] = c[0] - b[0];
	B[1] = c[1] - b[1];

	double dotAB = A[0] * B[0] + A[1] * B[1];
	double lengthA = sqrt(POW2(A[0]) + POW2(A[1]));
	double lengthB = sqrt(POW2(B[0]) + POW2(B[1]));
	double arg = dotAB / (lengthA * lengthB);
	arg = MAX(arg, -1.0);
	arg = MIN(arg, 1.0);

       	double angle;
	if (nb_utils2D_orient(a, b, c) > 0.0)
		angle = acos(arg);
	else
		angle = 2.0 * NB_PI - acos(arg);
	return angle;
}

void nb_utils2D_get_enveloping_box(uint32_t N_vertices,
				    const void* vertices,
				    uint16_t size_of_vtx_type,
				    double (*get_x)(const void *const vtx),
				    double (*get_y)(const void *const vtx),
				    double box[4])
{
	void *vtx = nb_array_get(vertices, size_of_vtx_type, 0);
	double x = get_x(vtx);
	double y = get_y(vtx);
	box[0] = x;
	box[1] = y;
	box[2] = x;
	box[3] = y;
	for (uint32_t i = 1; i < N_vertices; i++) {
		vtx = nb_array_get(vertices, size_of_vtx_type, i);
		x = get_x(vtx);
		y = get_y(vtx);
		if (x < box[0])
			box[0] = x;
		else if (x > box[2])
			box[2] = x;

		if (y < box[1])
			box[1] = y;
		else if (y > box[3])
			box[3] = y;
	}
}

void nb_utils2D_get_enveloping_box_from_subset
				(uint32_t N_subset,
				 const uint32_t subset[],
				 const void* vertices,
				 uint16_t size_of_vtx_type,
				 double (*get_x)(const void *const vtx),
				 double (*get_y)(const void *const vtx),
				 double box[4])
{	
	uint32_t id = subset[0];
	void *vtx = nb_array_get(vertices, size_of_vtx_type, id);
	double x = get_x(vtx);
	double y = get_y(vtx);
	box[0] = x;
	box[1] = y;
	box[2] = x;
	box[3] = y;
	for (uint32_t i = 1; i < N_subset; i++) {
		uint32_t id = subset[i];
		vtx = nb_array_get(vertices, size_of_vtx_type, id);
		x = get_x(vtx);
		y = get_y(vtx);
		if (x < box[0])
			box[0] = x;
		else if (x > box[2])
			box[2] = x;

		if (y < box[1])
			box[1] = y;
		else if (y > box[3])
			box[3] = y;
	}
}

double nb_utils2D_orient(const double t1[2],
				 const double t2[2],
				 const double t3[2])
{
	init_robust_predicates_if_not();
	return orient2d(t1, t2, t3);/* Robust predicate */
}

static void init_robust_predicates_if_not(void)
{
	if (init_id != 1) {
		exactinit();/* Init robust predicates */
		init_id = 1;
	}
}

bool nb_utils2D_is_in_half_side(const double v1[2],
				 const double v2[2],
				 const double v3[2])
{
	double sign = nb_utils2D_orient(v1, v2, v3);
	return (sign >= NB_GEOMETRIC_TOL);
}

double nb_utils2D_get_poly_area(const double *p[2], uint16_t N)
{	
	double orient = 0.0;
	for (uint16_t i = 0; i < N; i++) {
		double x1 = p[i][0];
		double y1 = p[i][1];
		uint16_t j = (i + 1) % N;
		double x2 = p[j][0];
		double y2 = p[j][1];
		orient += x1 * y2 - x2 * y1;
	}
	return 0.5 * orient;
}

double nb_utils2D_get_trg_area(const double t1[2],
				const double t2[2],
				const double t3[2])
{
	return 0.5 * nb_utils2D_orient(t1, t2, t3);
}

void nb_utils2D_trg_get_centroid(const double t1[2],
				  const double t2[2],
				  const double t3[2],
				  double centroid[2])
{
	centroid[0] = (t1[0] + t2[0] + t3[0]) / 3.0;
	centroid[1] = (t1[1] + t2[1] + t3[1]) / 3.0;
}

double nb_utils2D_get_circumradius(const double t1[2],
			const double t2[2],
			const double t3[2])
{
	const double Sk = nb_utils2D_orient(t1, t2, t3);
	/* Compute distance between vertices */
	const double L1 = nb_utils2D_get_dist(t1, t2);
	const double L2 = nb_utils2D_get_dist(t3, t2);
	const double L3 = nb_utils2D_get_dist(t3, t1);
	return (L1*L2*L3)/(2.0*Sk);
}

void nb_utils2D_get_circumcenter(const double t1[2],
				  const double t2[2],
				  const double t3[2],
				  double circumcenter[2])
{
	const double a = t3[0]-t1[0];
	const double b = t3[1]-t1[1];
	const double c = t3[0]-t2[0];
	const double d = t3[1]-t2[1];
	const double detA = a*d-b*c;
	const double ai = d/detA;
	const double bi = -b/detA;
	const double ci = -c/detA;
	const double di = a/detA;
	const double aux1 = POW2(t1[0]) + POW2(t1[1]);
	const double aux2 = POW2(t2[0]) + POW2(t2[1]);
	const double aux3 = POW2(t3[0]) + POW2(t3[1]);
	const double e = (aux3-aux1)/2.0;
	const double f = (aux3-aux2)/2.0;
	circumcenter[0] = ai*e + bi*f;
	circumcenter[1] = ci*e + di*f;
}

void nb_utils2D_get_circumcenter_from_sgm(const double s1[2],
					   const double s2[2],
					   double radii,
					   double circumcenter[2])
{
	double n[2];
	n[0] = s1[1] - s2[1];
	n[1] = s2[0] - s1[0];
	double dist = sqrt(POW2(n[0]) + POW2(n[1]));
	n[0] /= dist;
	n[1] /= dist;
	double a = sqrt(POW2(radii) - 0.25 * POW2(dist));
	circumcenter[0] = 0.5 * (s1[0] + s2[0]) + a * n[0];
	circumcenter[1] = 0.5 * (s1[1] + s2[1]) + a * n[1];
}

double nb_utils2D_get_trg_min_angle(const double t1[2],
				     const double t2[2],
				     const double t3[2])
{
	double cr2se = nb_utils2D_get_cr2se_ratio(t1, t2, t3);
	return asin(1.0 / (2.0 * cr2se));
}

double nb_utils2D_get_cr2se_ratio(const double t1[2],
				   const double t2[2],
				   const double t3[2])
/* Get circumradius to shortest edge ratio */
{
	double shortest_edge = nb_utils2D_get_min_trg_edge(t1, t2, t3);
	return nb_utils2D_get_circumradius(t1, t2, t3) / shortest_edge;
}

double nb_utils2D_get_trg_quality(const double t1[2],
				   const double t2[2],
				   const double t3[2])
{
	double B = nb_utils2D_get_cr2se_ratio(t1, t2, t3);
	return NB_MATH_INV_SQRT3/B;
}

double nb_utils2D_get_min_trg_edge(const double t1[2],
				    const double t2[2],
				    const double t3[2])
{
	const double l1 = nb_utils2D_get_dist2(t1, t2);
	const double l2 = nb_utils2D_get_dist2(t3, t2);
	const double l3 = nb_utils2D_get_dist2(t3, t1);
	double min = MIN(l1, l2);
	min = MIN(l3, min);
	return sqrt(min);
}

double nb_utils2D_get_max_trg_edge(const double t1[2],
				    const double t2[2],
				    const double t3[2])
{
	const double l1 = nb_utils2D_get_dist2(t1, t2);
	const double l2 = nb_utils2D_get_dist2(t3, t2);
	const double l3 = nb_utils2D_get_dist2(t3, t1);
	double max = MAX(l1, l2);
	max = MAX(l3, max);
	return sqrt(max);
}

uint32_t nb_utils2D_get_closest_vtx(double x, double y,
				     uint32_t N_vertices,
				     const double vertices[])
{
	double min_r = 2e30;
	uint32_t min_i = N_vertices + 1;
	for (uint32_t i = 0; i < N_vertices; i++) {
		double r = POW2(vertices[i*2] - x) +
			POW2(vertices[i*2+1] - y);
		if (r < min_r) {
			min_r = r;
			min_i = i;
		}
	}
	return min_i;
}

uint32_t nb_utils2D_get_closest_vtx_not_ignored
				(double x, double y,
				 uint32_t N_vertices,
				 const double vertices[],
				 uint32_t N_index_to_ignore,
				 const uint32_t index_to_ignore[])
{
	double min_r = 2e30;
	uint32_t min_i = N_vertices + 1;
	/* TEMPORAL: Use a mask array to make it faster */
	/* Slower version (ignoring some elements) */
	for (uint32_t i = 0; i < N_vertices; i++) {
		bool ignore = false;
		for (uint32_t j = 0; j  < N_index_to_ignore; j++) {
			if (i == index_to_ignore[j]) {
				ignore = true;
				break;
			}
		}
		if (ignore)
			continue;
		
		double r = POW2(vertices[i*2] - x) +
			POW2(vertices[i*2+1] - y);
		if (r < min_r) {
			min_r = r;
			min_i = i;
		}
	}
	return min_i;
}

void nb_utils2D_get_closest_pnt_to_sgm(const double s1[2],
					const double s2[2],
					const double p[2],
					double closest_point[2])
{
	const double u =
		((p[0]-s1[0])*(s2[0]-s1[0]) + (p[1]-s1[1])*(s2[1]-s1[1]))/
		((s2[0]-s1[0])*(s2[0]-s1[0]) + (s2[1]-s1[1])*(s2[1]-s1[1]));
	if (0 >= u) {
		memcpy(closest_point, s1, 2 * sizeof(*s1));
	} else if (1 <= u) {
		memcpy(closest_point, s2, 2 * sizeof(*s2));
	} else {
		closest_point[0] = s1[0] + u*(s2[0]-s1[0]);
		closest_point[1] = s1[1] + u*(s2[1]-s1[1]);
	}
}

nb_intersect_t nb_utils2D_get_sgm_intersection
				(const double a1[2], const double a2[2],
				 const double b1[2], const double b2[2],
				 /* Output NULL if not required */
				 double intersection[2])
{
	nb_intersect_t status;
	const double denominator = 
		((b2[1]-b1[1])*(a2[0]-a1[0]) - (b2[0]-b1[0])*(a2[1]-a1[1]));
	if (fabs(denominator) < NB_GEOMETRIC_TOL) {
		status = NB_PARALLEL;
		goto EXIT;
	}
	const double ua = ((b2[0]-b1[0])*(a1[1]-b1[1]) - 
			   (b2[1]-b1[1])*(a1[0]-b1[0])) /
		denominator;
	const double ub = ((a2[0]-a1[0])*(a1[1]-b1[1]) -
			   (a2[1]-a1[1])*(a1[0]-b1[0])) /
		denominator;

	if (ua < - NB_GEOMETRIC_TOL || ua > 1 + NB_GEOMETRIC_TOL) {
		status = NB_NOT_INTERSECTED;
		goto EXIT;
	}

	if (ub < - NB_GEOMETRIC_TOL || ub > 1 + NB_GEOMETRIC_TOL) {
		status = NB_NOT_INTERSECTED;
		goto EXIT;
	}

	/* Check extremes (for floating error) */
	double p[2];
	p[0] = a1[0] + ua*(a2[0]-a1[0]);
	p[1] = a1[1] + ua*(a2[1]-a1[1]);
	if (NULL != intersection) {
		intersection[0] = p[0];
		intersection[1] = p[1];
	}

	if (fabs(p[0]-a1[0]) < NB_GEOMETRIC_TOL &&
	    fabs(p[1]-a1[1]) < NB_GEOMETRIC_TOL) {
		status = NB_INTERSECT_ON_A1;
		goto EXIT;
	}
	if (fabs(p[0]-a2[0]) < NB_GEOMETRIC_TOL &&
	    fabs(p[1]-a2[1]) < NB_GEOMETRIC_TOL) {
		status = NB_INTERSECT_ON_A2;
		goto EXIT;
	}
	if (fabs(p[0]-b1[0]) < NB_GEOMETRIC_TOL &&
	    fabs(p[1]-b1[1]) < NB_GEOMETRIC_TOL) {
		status = NB_INTERSECT_ON_B1;
		goto EXIT;
	}
	if (fabs(p[0]-b2[0]) < NB_GEOMETRIC_TOL &&
	    fabs(p[1]-b2[1]) < NB_GEOMETRIC_TOL) {
		status = NB_INTERSECT_ON_B2;
		goto EXIT;
	}
	status = NB_INTERSECTED;
EXIT:
	return status;
}

bool nb_utils2D_are_sgm_intersected(const double a1[2], const double a2[2],
				     const double b1[2], const double b2[2],
				     /* Output NULL if not required */
				     double intersection[2])
{
	nb_intersect_t status =
		nb_utils2D_get_sgm_intersection(a1, a2, b1, b2,
						 intersection);
	return (status == NB_INTERSECTED);
}

bool nb_utils2D_sgm_intersects_trg(const double t1[2], const double t2[2],
				    const double t3[2], const double s1[2],
				    const double s2[2])
{
	nb_intersect_t status =
		nb_utils2D_get_sgm_intersection(t1, t2, s1, s2, NULL);
	if (NB_INTERSECTED != status) {
		status = nb_utils2D_get_sgm_intersection(t2, t3, s1, s2, NULL);
		if (NB_INTERSECTED != status) {
			status = nb_utils2D_get_sgm_intersection(t3, t1, 
								  s1, s2, NULL);
		}
	}
	return (NB_INTERSECTED == status);
}

bool nb_utils2D_sgm_intersects_circle(const double circumcenter[2],
				       double radius, const double s1[2],
				       const double s2[2])
{
	double p[2];
	nb_utils2D_get_closest_pnt_to_sgm(s1, s2, circumcenter, p);
	return (POW2(radius) - nb_utils2D_get_dist2(circumcenter, p) >
		NB_GEOMETRIC_TOL);
}

bool nb_utils2D_level_set_intersects_trg(double v1, double v2, double v3,
					 double level_set)
{
	bool out = true;
	if (level_set - v1 > -1e-6) {
		if (level_set - v2 > -1e-6)
			if (level_set - v3 > -1e-6)
				out = false;
	}
	if (out) {
		if (level_set - v1 < 1e-6)
			if (level_set - v2 < 1e-6)
				if (level_set - v3 < 1e-6)
					out = false;
	}
	return out;
}

void nb_utils2D_get_trg_level_set_intersection(const double t1[2],
					       const double t2[2],
					       const double t3[2],
					       double v1, double v2, double v3,
					       double level_set,
					       double a[2], double b[2])
{
	if (level_set_intersects_sgm(v1, v2, level_set)) {
		get_sgm_level_set_intersection(t1, t2, v1, v2, level_set, a);
		if (level_set_intersects_sgm(v2, v3, level_set))
			get_sgm_level_set_intersection(t2, t3, v2, v3,
						       level_set, b);
		else
			get_sgm_level_set_intersection(t1, t3, v1, v3,
						       level_set, b);
	} else {
		get_sgm_level_set_intersection(t2, t3, v2, v3, level_set, a);
		get_sgm_level_set_intersection(t1, t3, v1, v3, level_set, b);
	}
}

static inline bool level_set_intersects_sgm(double v1, double v2,
					    double level_set)
{
	return (level_set > v1 && level_set < v2) ||
		(level_set > v2 && level_set < v1);
}

static void get_sgm_level_set_intersection(const double s1[2],
					   const double s2[2],
					   double v1, double v2,
					   double level_set,
					   double p[2])
{
	double w = (level_set - v1) / (v2 - v1);
	p[0] = (1.0 - w) * s1[0] + w * s2[0];
	p[1] = (1.0 - w) * s1[1] + w * s2[1];
		
}

bool nb_utils2D_pnt_lies_on_sgm(const double s1[2], const double s2[2],
				 const double p[2])
{  
	if (fabs(s1[0] - s2[0]) > NB_GEOMETRIC_TOL &&
	    fabs(s1[1] - s2[1]) > NB_GEOMETRIC_TOL) {
		double wx = (p[0] - s2[0]) / (s1[0] - s2[0]);
		double wy = (p[1] - s2[1]) / (s1[1] - s2[1]);

		if (wx > 1.0 - NB_GEOMETRIC_TOL)
			return false;
		if (wx < NB_GEOMETRIC_TOL) 
			return false;
		if (fabs(wx - wy) > NB_GEOMETRIC_TOL)
			return false;
	} else if (fabs(s1[0] - s2[0]) > NB_GEOMETRIC_TOL) {
		double wx = (p[0] - s2[0]) / (s1[0] - s2[0]);      
		if (wx > 1.0 - NB_GEOMETRIC_TOL)
			return false;
		if (wx < NB_GEOMETRIC_TOL)
			return false;
		if (fabs(p[1] - s2[1]) > NB_GEOMETRIC_TOL)
			return false;
	} else {
		double wy = (p[1] - s2[1]) / (s1[1] - s2[1]);      
		if (wy > 1.0 - NB_GEOMETRIC_TOL)
			return false;
		if (wy < NB_GEOMETRIC_TOL) 
			return false;
		if (fabs(p[0] - s2[0]) > NB_GEOMETRIC_TOL)
			return false;
	}
	return true;
}

bool nb_utils2D_pnt_lies_in_trg(const double t1[2],
				 const double t2[2],
				 const double t3[2],
				 const double p [2])
{
	const double edge1 = nb_utils2D_orient(t1, t2, p);
	const double edge2 = nb_utils2D_orient(t2, t3, p);
	const double edge3 = nb_utils2D_orient(t3, t1, p);

	return (edge1 > -NB_GEOMETRIC_TOL &&
		edge2 > -NB_GEOMETRIC_TOL &&
		edge3 > -NB_GEOMETRIC_TOL);
}

bool nb_utils2D_pnt_lies_in_poly(int N, const double *poly,
				 const double p[2])
/* Edge crossing algorithm */
{
	int N_cross = 0;

	for (int i = 0; i < N; i++) {
		int next = (i+1) % N;
		if (right_ray_is_intersected(&(poly[i*2]),
					     &(poly[next*2]), p))
			N_cross += 1;
	}
	return N_cross & 1;
}

static bool right_ray_is_intersected(const double v1[2], const double v2[2],
				     const double p[2])
{
	bool out = true;
	if (v1[1] > p[1] && v2[1] > p[1]) {
		out = false;
		goto EXIT;
	}
	if (v1[1] < p[1] && v2[1] < p[1]) {
		out = false;
		goto EXIT;
	}
	if (v1[0] < p[0] && v2[0] < p[0]) {
		out = false;
		goto EXIT;
	}
	if (v1[0] + v2[0] - 2 * p[0] < 0) {
		out = false;
		goto EXIT;
	}

EXIT:
	return out;
}


bool nb_utils2D_pnt_lies_in_poly_bnd(int N, const double *poly,
				     const double p[2])
/* Edge crossing algorithm */
{
	bool in = false;
	for (int i = 0; i < N; i++) {
		int next = (i+1) % N;
		in = nb_utils2D_pnt_lies_on_sgm(&(poly[i*2]),
						&(poly[next*2]), p);
		if (in)
			break;
	}
	return in;
}

bool nb_utils2D_pnt_lies_in_diametral_circle(const double s1[2],
					      const double s2[2],
					      const double p[2])
{
	return NB_GEOMETRIC_TOL >
		(s1[0] - p[0]) * (s2[0] - p[0]) +
		(s1[1] - p[1]) * (s2[1] - p[1]);
}

bool nb_utils2D_pnt_lies_strictly_in_diametral_circle(const double s1[2],
						       const double s2[2],
						       const double p[2])
{
	return -NB_GEOMETRIC_TOL >
		(s1[0] - p[0]) * (s2[0] - p[0]) +
		(s1[1] - p[1]) * (s2[1] - p[1]);
}

bool nb_utils2D_pnt_lies_strictly_in_circumcircle(const double t1[2],
						   const double t2[2],
						   const double t3[2],
						   const double p[2])
{
	double det = det_circumcircle(t1, t2, t3, p);
	return (det > INCIRCLE_TOLERANCE);
}

static double det_circumcircle(const double t1[2],
			       const double t2[2],
			       const double t3[2],
			       const double p[2])
{
	init_robust_predicates_if_not();
	return incircle(t1, t2, t3, p);/* Robust predicate */
}

bool nb_utils2D_pnt_is_cocircular(const double t1[2],
				  const double t2[2],
				  const double t3[2],
				  const double p[2])
{
	double det = det_circumcircle(t1, t2, t3, p);
	return (fabs(det) < INCIRCLE_TOLERANCE);
}

bool nb_utils2D_pnt_lies_in_box(const double box[4],
				 const double p[2])
{
	if (p[0] < box[0])
		return false;
	if (p[0] > box[2])
		return false;
	if (p[1] < box[1])
		return false;
	if (p[1] > box[3])
		return false;
	return true;
}
