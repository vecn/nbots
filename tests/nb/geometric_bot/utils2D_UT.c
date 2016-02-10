#include <stdbool.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "nb/geometric_bot/utils2D.h"

#include "test_library.h"
#include "test_add.h"

#define TOLERANCE (1e-9)
#define INV_SQRT2 (0.707106781187)
#define INV_SQRT3 (0.577350269190)
#define SQRT2 (1.41421356237)
#define PI (3.14159265359)

#define POW2(a) ((a)*(a))

static bool check_get_x_from_darray(void);
static bool check_get_y_from_darray(void);
static bool check_get_normal(void);
static bool check_get_dist(void);
static bool check_get_dist2(void);
static bool check_get_delaunay_dist(void);
static bool check_get_delaunay_dist_neg(void);
static bool check_get_enveloping_box(void);
static bool check_get_enveloping_box_from_subset(void);
static bool check_get_2x_trg_area(void);
static bool check_get_trg_area(void);
static bool check_get_trg_centroid(void);
static bool check_get_circumradius(void);
static bool check_get_circumcenter(void);
static bool check_get_circumcenter_from_sgm(void);
static bool check_get_trg_min_angle(void);
static bool check_get_cr2se_ratio(void);
static bool check_get_trg_quality(void);
static bool check_get_min_trg_edge(void);
static bool check_get_max_trg_edge(void);
static bool check_get_closest_vtx(void);
static bool check_get_closest_vtx_not_ignored(void);
static bool check_get_closest_pnt_to_sgm_case_middle(void);
static bool check_get_closest_pnt_to_sgm_case_vtx1(void);
static bool check_get_closest_pnt_to_sgm_case_vtx2(void);
static bool check_are_sgm_intersected_true(void);
static bool check_are_sgm_intersected_parallel(void);
static bool check_are_sgm_intersected_strictly_false(void);
static bool check_are_sgm_intersected_a1(void);
static bool check_are_sgm_intersected_a2(void);
static bool check_are_sgm_intersected_b1(void);
static bool check_are_sgm_intersected_b2(void);
static bool check_sgm_intersects_trg_1edge(void);
static bool check_sgm_intersects_trg_2edge(void);
static bool check_sgm_intersects_trg_false(void);
static bool check_sgm_intersects_circle_1x(void);
static bool check_sgm_intersects_circle_2x(void);
static bool check_sgm_intersects_circle_false(void);
static bool check_pnt_lies_on_sgm_true(void);
static bool check_pnt_lies_on_sgm_false(void);
static bool check_pnt_lies_in_trg_true(void);
static bool check_pnt_lies_in_trg_true_on_limit(void);
static bool check_pnt_lies_in_trg_false(void);
static bool check_pnt_lies_strictly_in_trg_true(void);
static bool check_pnt_lies_strictly_in_trg_false_on_limit(void);
static bool check_pnt_lies_strictly_in_trg_false(void);
static bool check_pnt_lies_strictly_in_diametral_circle_true(void);
static bool check_pnt_lies_strictly_in_diametral_circle_false_on_limit(void);
static bool check_pnt_lies_strictly_in_diametral_circle_false(void);
static bool check_pnt_lies_strictly_in_circumcircle_true(void);
static bool check_pnt_lies_strictly_in_circumcircle_false_on_limit(void);
static bool check_pnt_lies_strictly_in_circumcircle_false(void);
static bool check_pnt_lies_in_box_true(void);
static bool check_pnt_lies_in_box_true_on_limit(void);
static bool check_pnt_lies_in_box_false(void);

inline int vcn_test_get_driver_id(void)
{
	return NB_DRIVER_UNIT_TEST;
}

void vcn_test_load_tests(void *tests_ptr)
{
	vcn_test_add(tests_ptr, check_get_x_from_darray,
		     "Check get_x_from_darray()");
	vcn_test_add(tests_ptr, check_get_y_from_darray,
		     "Check get_y_from_darray()");
	vcn_test_add(tests_ptr, check_get_normal,
		     "Check get_normal()");
	vcn_test_add(tests_ptr, check_get_dist,
		     "Check get_dist()");
	vcn_test_add(tests_ptr, check_get_dist2,
		     "Check get_dist2()");
	vcn_test_add(tests_ptr, check_get_delaunay_dist,
		     "Check get_delaunay_dist()");
	vcn_test_add(tests_ptr, check_get_delaunay_dist_neg,
		     "Check get_delaunay_dist_neg()");
	vcn_test_add(tests_ptr, check_get_enveloping_box,
		     "Check get_enveloping_box()");
	vcn_test_add(tests_ptr, check_get_enveloping_box_from_subset,
		     "Check get_enveloping_box_from_subset()");
	vcn_test_add(tests_ptr, check_get_2x_trg_area,
		     "Check get_2x_trg_area()");
	vcn_test_add(tests_ptr, check_get_trg_area,
		     "Check get_trg_area()");
	vcn_test_add(tests_ptr, check_get_trg_centroid,
		     "Check get_trg_centroid()");
	vcn_test_add(tests_ptr, check_get_circumradius,
		     "Check get_circumradius()");
	vcn_test_add(tests_ptr, check_get_circumcenter,
		     "Check get_circumcenter()");
	vcn_test_add(tests_ptr, check_get_circumcenter_from_sgm,
		     "Check get_circumcenter_from_sgm()");
	vcn_test_add(tests_ptr, check_get_trg_min_angle,
		     "Check get_trg_min_angle()");
	vcn_test_add(tests_ptr, check_get_cr2se_ratio,
		     "Check get_cr2se_ratio()");
	vcn_test_add(tests_ptr, check_get_trg_quality,
		     "Check get_trg_quality()");
	vcn_test_add(tests_ptr, check_get_min_trg_edge,
		     "Check get_min_trg_edge()");
	vcn_test_add(tests_ptr, check_get_max_trg_edge,
		     "Check get_max_trg_edge()");
	vcn_test_add(tests_ptr, check_get_closest_vtx,
		     "Check get_closest_vtx()");
	vcn_test_add(tests_ptr, check_get_closest_vtx_not_ignored,
		     "Check get_closest_vtx_not_ignored()");
	vcn_test_add(tests_ptr, check_get_closest_pnt_to_sgm_case_middle,
		     "Check get_closest_pnt_to_sgm() if closest is in sgm");
	vcn_test_add(tests_ptr, check_get_closest_pnt_to_sgm_case_vtx1,
		     "Check get_closest_pnt_to_sgm() if closest is vtx 1");
	vcn_test_add(tests_ptr, check_get_closest_pnt_to_sgm_case_vtx2,
		     "Check get_closest_pnt_to_sgm() if closest is vtx 2");
	vcn_test_add(tests_ptr, check_are_sgm_intersected_true,
		     "Check are_sgm_intersected() with true intersection");
	vcn_test_add(tests_ptr, check_are_sgm_intersected_parallel,
		     "Check are_sgm_intersected() if sgm are parallel");
	vcn_test_add(tests_ptr, check_are_sgm_intersected_strictly_false,
		     "Check are_sgm_intersected() if sgm are not intersected");
	vcn_test_add(tests_ptr, check_are_sgm_intersected_a1,
		     "Check are_sgm_intersected() if a1 is the intersection");
	vcn_test_add(tests_ptr, check_are_sgm_intersected_a2,
		     "Check are_sgm_intersected() if a2 is the intersection");
	vcn_test_add(tests_ptr, check_are_sgm_intersected_b1,
		     "Check are_sgm_intersected() if b1 is the intersection");
	vcn_test_add(tests_ptr, check_are_sgm_intersected_b2,
		     "Check are_sgm_intersected() if b2 is the intersection");
	vcn_test_add(tests_ptr, check_sgm_intersects_trg_1edge,
		     "Check sgm_intersects_trg() with one edge intersection");
	vcn_test_add(tests_ptr, check_sgm_intersects_trg_2edge,
		     "Check sgm_intersects_trg() with two edges intersection");
	vcn_test_add(tests_ptr, check_sgm_intersects_trg_false,
		     "Check sgm_intersects_trg() with no intersections");
	vcn_test_add(tests_ptr, check_sgm_intersects_circle_1x,
		     "Check sgm_intersects_circle() with one intersection");
	vcn_test_add(tests_ptr, check_sgm_intersects_circle_2x,
		     "Check sgm_intersects_circle() with two intersection");
	vcn_test_add(tests_ptr, check_sgm_intersects_circle_false,
		     "Check sgm_intersects_circle() with no intersections");
	vcn_test_add(tests_ptr, check_pnt_lies_on_sgm_true,
		     "Check pnt_lies_on_sgm() if lies inside");
	vcn_test_add(tests_ptr, check_pnt_lies_on_sgm_false,
		     "Check pnt_lies_on_sgm() if lies outside");
	vcn_test_add(tests_ptr, check_pnt_lies_in_trg_true,
		     "Check pnt_lies_in_trg() if lies inside");
	vcn_test_add(tests_ptr, check_pnt_lies_in_trg_true_on_limit,
		     "Check pnt_lies_in_trg() if lies in the edge");
	vcn_test_add(tests_ptr, check_pnt_lies_in_trg_false,
		     "Check pnt_lies_in_trg() if lies outside");
	vcn_test_add(tests_ptr, check_pnt_lies_strictly_in_trg_true,
		     "Check pnt_lies_strictly_in_trg() if lies inside");
	vcn_test_add(tests_ptr, check_pnt_lies_strictly_in_trg_false_on_limit,
		     "Check pnt_lies_strictly_in_trg() if lies in the edge");
	vcn_test_add(tests_ptr, check_pnt_lies_strictly_in_trg_false,
		     "Check pnt_lies_strictly_in_trg() if lies outside");
	vcn_test_add(tests_ptr,
		     check_pnt_lies_strictly_in_diametral_circle_true,
		     "Check pnt_lies_strict_in_diam_circ() if lies inside");
	vcn_test_add(tests_ptr,
		     check_pnt_lies_strictly_in_diametral_circle_false_on_limit,
		     "Check pnt_lies_strict_in_diam_circ() if lies on limits");
	vcn_test_add(tests_ptr,
		     check_pnt_lies_strictly_in_diametral_circle_false,
		     "Check pnt_lies_strict_in_diam_circ() if lies outside");
	vcn_test_add(tests_ptr, check_pnt_lies_strictly_in_circumcircle_true,
		     "Check pnt_lies_strict_in_circumc() if lies inside");
	vcn_test_add(tests_ptr,
		     check_pnt_lies_strictly_in_circumcircle_false_on_limit,
		     "Check pnt_lies_strict_in_circumc() if lies on limits");
	vcn_test_add(tests_ptr, check_pnt_lies_strictly_in_circumcircle_false,
		     "Check pnt_lies_strict_in_circumc() if lies outside");
	vcn_test_add(tests_ptr, check_pnt_lies_in_box_true,
		     "Check pnt_lies_in_box() if lies inside");
	vcn_test_add(tests_ptr, check_pnt_lies_in_box_true_on_limit,
		     "Check pnt_lies_in_box() if lies on limits");
	vcn_test_add(tests_ptr, check_pnt_lies_in_box_false,
		     "Check pnt_lies_in_box() if lies outside");
}

static bool check_get_x_from_darray(void)
{
	int id = 2;
	double vertices[10] = {0, 5, 2, 1, 8, 9, 6, 3, 7, 4};
	void *vtx = &(vertices[id * 2]);
	double x = vcn_utils2D_get_x_from_darray(vtx);
	return (vertices[id * 2] == x);
}

static bool check_get_y_from_darray(void)
{
	int id = 3;
	double vertices[10] = {0, 5, 2, 1, 8, 9, 6, 3, 7, 4};
	void *vtx = &(vertices[id * 2]);
	double y = vcn_utils2D_get_y_from_darray(vtx);
	return (vertices[id*2+1] == y);
}

static bool check_get_normal(void)
{
	double x1[2] = {0, 0};
	double x2[2] = {1, 1};
	double normal[2];
	vcn_utils2D_get_normal(x1, x2, normal);
	double dist = sqrt(POW2(normal[0]) + POW2(normal[1]));
	return (fabs(1.0 - dist) < TOLERANCE) &&
		(fabs(normal[0] - INV_SQRT2) < TOLERANCE) &&
		(fabs(normal[1] - INV_SQRT2) < TOLERANCE);
}

static bool check_get_dist(void)
{
	double p1[2] = {0, 0};
	double p2[2] = {1, 1};
	double dist = vcn_utils2D_get_dist(p1, p2);
	return fabs(dist - SQRT2) < TOLERANCE;
}

static bool check_get_dist2(void)
{
	double p1[2] = {0, 0};
	double p2[2] = {1, 1};
	double dist = vcn_utils2D_get_dist2(p1, p2);
	return fabs(dist - 2) < TOLERANCE;
}

static bool check_get_delaunay_dist(void)
{
	double p1[2] = {0, 0};
	double p2[2] = {1, 0};
	double p3[2] = {1, 1};
	double dist = vcn_utils2D_get_delaunay_dist(p1, p2, p3);
	return fabs(dist - 0.5) < TOLERANCE;
}

static bool check_get_delaunay_dist_neg(void)
{
	double p1[2] = {0, 0};
	double p2[2] = {2, 0};
	double p3[2] = {1, 0.5};
	double dist = vcn_utils2D_get_delaunay_dist(p1, p2, p3);
	return dist < 0.0;
}

static bool check_get_enveloping_box(void)
{
	double vertices[10] = {0, 0, 0, 1, 1, 1, 1, 0, 0.5, 0.5};
	int16_t vtx_size = 2 * sizeof(double);
	double box[4];
	vcn_utils2D_get_enveloping_box(5, vertices, vtx_size,
				       vcn_utils2D_get_x_from_darray,
				       vcn_utils2D_get_y_from_darray,
				       box);
	return (fabs(box[0] - 0) < TOLERANCE) &&
		(fabs(box[1] - 0) < TOLERANCE) &&
		(fabs(box[2] - 1) < TOLERANCE) &&
		(fabs(box[3] - 1) < TOLERANCE);
}

static bool check_get_enveloping_box_from_subset(void)
{
	double vertices[10] = {0, 0, 0, 1, 1, 1, 1, 0, 0.5, 0.5};
	int16_t vtx_size = 2 * sizeof(double);
	uint32_t subset[2] = {0, 4};
	double box[4];
	vcn_utils2D_get_enveloping_box_from_subset(2, subset,
						   vertices, vtx_size,
						   vcn_utils2D_get_x_from_darray,
						   vcn_utils2D_get_y_from_darray,
						   box);
	return (fabs(box[0] - 0) < TOLERANCE) &&
		(fabs(box[1] - 0) < TOLERANCE) &&
		(fabs(box[2] - 0.5) < TOLERANCE) &&
		(fabs(box[3] - 0.5) < TOLERANCE);
}

static bool check_get_2x_trg_area(void)
{
	double t1[2] = {0, 0};
	double t2[2] = {1, 0};
	double t3[2] = {1, 1};
	double area_x2 = vcn_utils2D_get_2x_trg_area(t1, t2, t3);
	return fabs(area_x2 - 1) < TOLERANCE;
}

static bool check_get_trg_area(void)
{
	double t1[2] = {0, 0};
	double t2[2] = {1, 0};
	double t3[2] = {1, 1};
	double area = vcn_utils2D_get_trg_area(t1, t2, t3);
	return fabs(area - 0.5) < TOLERANCE;
}

static bool check_get_trg_centroid(void)
{
	double t1[2] = {0, 0};
	double t2[2] = {1, 0};
	double t3[2] = {1, 1};
	double centroid[2];
	vcn_utils2D_get_trg_centroid(t1, t2, t3, centroid);
	return (fabs(centroid[0] - 2.0/3.0) < TOLERANCE) &&
		(fabs(centroid[1] - 1.0/3.0) < TOLERANCE);
}

static bool check_get_circumradius(void)
{
	double t1[2] = {0, 0};
	double t2[2] = {1, 0};
	double t3[2] = {1, 1};
	double circumradius = vcn_utils2D_get_circumradius(t1, t2, t3);
	return (fabs(circumradius - SQRT2/2.0) < TOLERANCE);
}

static bool check_get_circumcenter(void)
{
	double t1[2] = {0, 0};
	double t2[2] = {1, 0};
	double t3[2] = {1, 1};
	double circumcenter[2];
	vcn_utils2D_get_circumcenter(t1, t2, t3, circumcenter);
	return (fabs(circumcenter[0] - 0.5) < TOLERANCE) &&
		(fabs(circumcenter[1] - 0.5) < TOLERANCE);
}

static bool check_get_circumcenter_from_sgm(void)
{
	double s1[2] = {0, 0};
	double s2[2] = {1, 0};
	double circumcenter[2];
	vcn_utils2D_get_circumcenter_from_sgm(s1, s2, SQRT2/2.0, circumcenter);
	return (fabs(circumcenter[0] - 0.5) < TOLERANCE) &&
		(fabs(circumcenter[1] - 0.5) < TOLERANCE);
}

static bool check_get_trg_min_angle(void)
{
	double t1[2] = {0, 0};
	double t2[2] = {1, 0};
	double t3[2] = {1, 1};
	double min_angle = vcn_utils2D_get_trg_min_angle(t1, t2, t3);
	return (fabs(min_angle - PI/4.0) < TOLERANCE);
}

static bool check_get_cr2se_ratio(void)
{
	double t1[2] = {0, 0};
	double t2[2] = {1, 0};
	double t3[2] = {1, 1};
	double cr2se = vcn_utils2D_get_cr2se_ratio(t1, t2, t3);
	return (fabs(cr2se - SQRT2/2) < TOLERANCE);
}

static bool check_get_trg_quality(void)
{
	double t1[2] = {0, 0};
	double t2[2] = {1, 0};
	double t3[2] = {1, 1};
	double quality = vcn_utils2D_get_trg_quality(t1, t2, t3);
	double cr2se = vcn_utils2D_get_cr2se_ratio(t1, t2, t3);
	return (fabs(quality - INV_SQRT3/cr2se) < TOLERANCE);
}

static bool check_get_min_trg_edge(void)
{
	double t1[2] = {0, 0};
	double t2[2] = {1, 0};
	double t3[2] = {1, 1};
	double min = vcn_utils2D_get_min_trg_edge(t1, t2, t3);
	return (fabs(min - 1.0) < TOLERANCE);
}

static bool check_get_max_trg_edge(void)
{
	double t1[2] = {0, 0};
	double t2[2] = {1, 0};
	double t3[2] = {1, 1};
	double max = vcn_utils2D_get_max_trg_edge(t1, t2, t3);
	return (fabs(max - SQRT2) < TOLERANCE);
}

static bool check_get_closest_vtx(void)
{
	double vertices[10] = {0, 0, 0, 1, 1, 1, 1, 0, 0.5, 0.5};
	uint32_t id_closest = vcn_utils2D_get_closest_vtx(2, 2,
							  5, vertices);
	return (2 == id_closest);

}

static bool check_get_closest_vtx_not_ignored(void)
{
	double vertices[10] = {0, 0, 0, 1, 1, 1, 1, 0, 0.5, 0.5};
	uint32_t ignore[2] = {0, 2};
	uint32_t id_closest =
		vcn_utils2D_get_closest_vtx_not_ignored(-1, -1,
							5, vertices,
							2, ignore);
	return (4 == id_closest);

}

static bool check_get_closest_pnt_to_sgm_case_middle(void)
{
	double s1[2] = {0, 0};
	double s2[2] = {1, 1};
	double p[2] = {0, 1};
	double cp[2];
	vcn_utils2D_get_closest_pnt_to_sgm(s1, s2, p, cp);
	return (fabs(cp[0] - 0.5) < TOLERANCE) &&
		(fabs(cp[1] - 0.5) < TOLERANCE);
}

static bool check_get_closest_pnt_to_sgm_case_vtx1(void)
{
	double s1[2] = {0, 0};
	double s2[2] = {1, 1};
	double p[2] = {0, -1};
	double cp[2];
	vcn_utils2D_get_closest_pnt_to_sgm(s1, s2, p, cp);
	return (fabs(cp[0] - 0) < TOLERANCE) &&
		(fabs(cp[1] - 0) < TOLERANCE);
}

static bool check_get_closest_pnt_to_sgm_case_vtx2(void)
{
	double s1[2] = {0, 0};
	double s2[2] = {1, 1};
	double p[2] = {1, 2};
	double cp[2];
	vcn_utils2D_get_closest_pnt_to_sgm(s1, s2, p, cp);
	return (fabs(cp[0] - 1) < TOLERANCE) &&
		(fabs(cp[1] - 1) < TOLERANCE);
}

static bool check_are_sgm_intersected_true(void)
{
	double a1[2] = {0, 0};
	double a2[2] = {1, 1};
	double b1[2] = {1, 0};
	double b2[2] = {0, 1};
	double intersection[2];
	int status;
	bool intersected = 
		vcn_utils2D_are_sgm_intersected(a1, a2, b1, b2,
						intersection, &status);
	return intersected &&
		(fabs(intersection[0] - 0.5) < TOLERANCE) &&
		(fabs(intersection[1] - 0.5) < TOLERANCE) &&
		(0 == status);
}

static bool check_are_sgm_intersected_parallel(void)
{
	double a1[2] = {0, 0};
	double a2[2] = {0, 1};
	double b1[2] = {1, 0};
	double b2[2] = {1, 1};
	double intersection[2];
	int status;
	bool intersected = 
		vcn_utils2D_are_sgm_intersected(a1, a2, b1, b2,
						intersection, &status);
	return !intersected && (3 == status);
}

static bool check_are_sgm_intersected_strictly_false(void)
{
	double a1[2] = {0, 0};
	double a2[2] = {1, 1};
	double b1[2] = {2, 1};
	double b2[2] = {1, 2};
	double intersection[2];
	int status;
	bool intersected = 
		vcn_utils2D_are_sgm_intersected(a1, a2, b1, b2,
						intersection, &status);
	return !intersected && (1 == status);
}

static bool check_are_sgm_intersected_a1(void)
{
	double a1[2] = {0, 0};
	double a2[2] = {1, 1};
	double b1[2] = {-1, 1};
	double b2[2] = {1, -1};
	double intersection[2];
	int status;
	bool intersected = 
		vcn_utils2D_are_sgm_intersected(a1, a2, b1, b2,
						intersection, &status);
	return !intersected &&
		(fabs(intersection[0] - a1[0]) < TOLERANCE) &&
		(fabs(intersection[1] - a1[1]) < TOLERANCE) &&
		(4 == status);
}

static bool check_are_sgm_intersected_a2(void)
{
	double a1[2] = {1, 1};
	double a2[2] = {0, 0};
	double b1[2] = {-1, 1};
	double b2[2] = {1, -1};
	double intersection[2];
	int status;
	bool intersected = 
		vcn_utils2D_are_sgm_intersected(a1, a2, b1, b2,
						intersection, &status);
	return !intersected &&
		(fabs(intersection[0] - a2[0]) < TOLERANCE) &&
		(fabs(intersection[1] - a2[1]) < TOLERANCE) &&
		(5 == status);
}

static bool check_are_sgm_intersected_b1(void)
{
	double a1[2] = {-1, 1};
	double a2[2] = {1, -1};
	double b1[2] = {0, 0};
	double b2[2] = {1, 1};
	double intersection[2];
	int status;
	bool intersected = 
		vcn_utils2D_are_sgm_intersected(a1, a2, b1, b2,
						intersection, &status);
	return !intersected &&
		(fabs(intersection[0] - b1[0]) < TOLERANCE) &&
		(fabs(intersection[1] - b1[1]) < TOLERANCE) &&
		(6 == status);
}

static bool check_are_sgm_intersected_b2(void)
{
	double a1[2] = {-1, 1};
	double a2[2] = {1, -1};
	double b1[2] = {1, 1};
	double b2[2] = {0, 0};
	double intersection[2];
	int status;
	bool intersected = 
		vcn_utils2D_are_sgm_intersected(a1, a2, b1, b2,
						intersection, &status);
	return !intersected &&
		(fabs(intersection[0] - b2[0]) < TOLERANCE) &&
		(fabs(intersection[1] - b2[1]) < TOLERANCE) &&
		(7 == status);
}

static bool check_sgm_intersects_trg_1edge(void)
{
	double t1[2] = {0, 0};
	double t2[2] = {1, 0};
	double t3[2] = {1, 1};
	double s1[2] = {1, 0};
	double s2[2] = {0.2, 0.8};
	bool intersects = vcn_utils2D_sgm_intersects_trg(t1, t2, t3, s1, s2);
	return intersects;
}

static bool check_sgm_intersects_trg_2edge(void)
{
	double t1[2] = {0, 0};
	double t2[2] = {1, 0};
	double t3[2] = {1, 1};
	double s1[2] = {1, 0};
	double s2[2] = {0.2, 1.2};
	bool intersects = vcn_utils2D_sgm_intersects_trg(t1, t2, t3, s1, s2);
	return intersects;
}

static bool check_sgm_intersects_trg_false(void)
{
	double t1[2] = {0, 0};
	double t2[2] = {1, 0};
	double t3[2] = {1, 1};
	double s1[2] = {1, 0};
	double s2[2] = {0.2, -0.8};
	bool not_intersects = 
		!vcn_utils2D_sgm_intersects_trg(t1, t2, t3, s1, s2);
	return not_intersects;
}

static bool check_sgm_intersects_circle_1x(void)
{
	double s1[2] = {0, 0};
	double s2[2] = {1, 1};
	double cc[2] = {1, 1};
	bool intersects = vcn_utils2D_sgm_intersects_circle(cc, 0.5, s1, s2);
	return intersects;
}

static bool check_sgm_intersects_circle_2x(void)
{
	double s1[2] = {0, 0};
	double s2[2] = {2, 2};
	double cc[2] = {1, 1};
	bool intersects = vcn_utils2D_sgm_intersects_circle(cc, 0.5, s1, s2);
	return intersects;
}

static bool check_sgm_intersects_circle_false(void)
{
	double s1[2] = {0, 0};
	double s2[2] = {-1, -1};
	double cc[2] = {1, 1};
	bool not_intersects = 
		!vcn_utils2D_sgm_intersects_circle(cc, 0.5, s1, s2);
	return not_intersects;
}

static bool check_pnt_lies_on_sgm_true(void)
{
	double s1[2] = {0, 0};
	double s2[2] = {1, 1};
	double p[2] = {0.5, 0.5};
	bool lies_inside = vcn_utils2D_pnt_lies_on_sgm(s1, s2, p);
	return lies_inside;
}

static bool check_pnt_lies_on_sgm_false(void)
{
	double s1[2] = {0, 0};
	double s2[2] = {1, 1};
	double p[2] = {1.1, 1.1};
	bool lies_outside = !vcn_utils2D_pnt_lies_on_sgm(s1, s2, p);
	return lies_outside;
}

static bool check_pnt_lies_in_trg_true(void)
{
	double t1[2] = {0, 0};
	double t2[2] = {1, 0};
	double t3[2] = {1, 1};
	double p[2] = {0.6, 0.5};
	bool lies_inside = vcn_utils2D_pnt_lies_in_trg(t1, t2, t3, p);
	return lies_inside;
}

static bool check_pnt_lies_in_trg_true_on_limit(void)
{
	double t1[2] = {0, 0};
	double t2[2] = {1, 0};
	double t3[2] = {1, 1};
	double p[2] = {0.5, 0.5};
	bool lies_inside = vcn_utils2D_pnt_lies_in_trg(t1, t2, t3, p);
	return lies_inside;
}

static bool check_pnt_lies_in_trg_false(void)
{
	double t1[2] = {0, 0};
	double t2[2] = {1, 0};
	double t3[2] = {1, 1};
	double p[2] = {0.4, 0.5};
	bool lies_outside = !vcn_utils2D_pnt_lies_in_trg(t1, t2, t3, p);
	return lies_outside;
}

static bool check_pnt_lies_strictly_in_trg_true(void)
{
	double t1[2] = {0, 0};
	double t2[2] = {1, 0};
	double t3[2] = {1, 1};
	double p[2] = {0.6, 0.5};
	bool lies_inside = vcn_utils2D_pnt_lies_strictly_in_trg(t1, t2, t3, p);
	return lies_inside;
}

static bool check_pnt_lies_strictly_in_trg_false_on_limit(void)
{
	double t1[2] = {0, 0};
	double t2[2] = {1, 0};
	double t3[2] = {1, 1};
	double p[2] = {0.5, 0.5};
	bool lies_outside =
		!vcn_utils2D_pnt_lies_strictly_in_trg(t1, t2, t3, p);
	return lies_outside;
}

static bool check_pnt_lies_strictly_in_trg_false(void)
{
	double t1[2] = {0, 0};
	double t2[2] = {1, 0};
	double t3[2] = {1, 1};
	double p[2] = {0.4, 0.5};
	bool lies_outside =
		!vcn_utils2D_pnt_lies_strictly_in_trg(t1, t2, t3, p);
	return lies_outside;
}

static bool check_pnt_lies_strictly_in_diametral_circle_true(void)
{
	double s1[2] = {0, 0};
	double s2[2] = {1, 0};
	double p[2] = {0.5, 0};
	bool lies_inside =
		vcn_utils2D_pnt_lies_strictly_in_diametral_circle(s1, s2, p);
	return lies_inside;
}

static bool check_pnt_lies_strictly_in_diametral_circle_false_on_limit(void)
{
	double s1[2] = {0, 0};
	double s2[2] = {1, 0};
	double p[2] = {0.5, 0.5};
	bool lies_outside =
		!vcn_utils2D_pnt_lies_strictly_in_diametral_circle(s1, s2, p);
	return lies_outside;
}

static bool check_pnt_lies_strictly_in_diametral_circle_false(void)
{
	double s1[2] = {0, 0};
	double s2[2] = {1, 0};
	double p[2] = {1, 0.5};
	bool lies_outside =
		!vcn_utils2D_pnt_lies_strictly_in_diametral_circle(s1, s2, p);
	return lies_outside;
}

static bool check_pnt_lies_strictly_in_circumcircle_true(void)
{
	double t1[2] = {0, 0};
	double t2[2] = {1, 0};
	double t3[2] = {1, 1};
	double p[2] = {0.5, 0.5};
	bool lies_inside =
		vcn_utils2D_pnt_lies_strictly_in_circumcircle(t1, t2, t3, p);
	return lies_inside;
}

static bool check_pnt_lies_strictly_in_circumcircle_false_on_limit(void)
{
	double t1[2] = {0, 0};
	double t2[2] = {1, 0};
	double t3[2] = {1, 1};
	double p[2] = {1, 0};
	bool lies_outside =
		!vcn_utils2D_pnt_lies_strictly_in_circumcircle(t1, t2, t3, p);
	return lies_outside;
}

static bool check_pnt_lies_strictly_in_circumcircle_false(void)
{
	double t1[2] = {0, 0};
	double t2[2] = {1, 0};
	double t3[2] = {1, 1};
	double p[2] = {2, 0};
	bool lies_outside =
		!vcn_utils2D_pnt_lies_strictly_in_circumcircle(t1, t2, t3, p);
	return lies_outside;
}

static bool check_pnt_lies_in_box_true(void)
{
	double box[4] = {0, 0, 1, 1};
	double p[2] = {0.5, 0.5};
	bool lies_inside = vcn_utils2D_pnt_lies_in_box(box, p);
	return lies_inside;
}

static bool check_pnt_lies_in_box_true_on_limit(void)
{
	double box[4] = {0, 0, 1, 1};
	double p[2] = {0.5, 1};
	bool lies_inside = vcn_utils2D_pnt_lies_in_box(box, p);
	return lies_inside;
}

static bool check_pnt_lies_in_box_false(void)
{
	double box[4] = {0, 0, 1, 1};
	double p[2] = {0.5, 1.1};
	bool lies_outside = !vcn_utils2D_pnt_lies_in_box(box, p);
	return lies_outside;
}
