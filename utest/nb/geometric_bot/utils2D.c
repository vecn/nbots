#include <stdbool.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include <CUnit/Basic.h>

#include "nb/geometric_bot/utils2D.h"

#define TOLERANCE (1e-9)
#define INV_SQRT2 (0.707106781187)
#define INV_SQRT3 (0.577350269190)
#define SQRT2 (1.41421356237)
#define PI (3.14159265359)

#define POW2(a) ((a)*(a))

static int suite_init(void);
static int suite_clean(void);

static void test_get_x_from_darray(void);
static void test_get_y_from_darray(void);
static void test_get_normal(void);
static void test_get_dist(void);
static void test_get_dist2(void);
static void test_get_delaunay_dist(void);
static void test_get_delaunay_dist_neg(void);
static void test_get_enveloping_box(void);
static void test_get_enveloping_box_from_subset(void);
static void test_orient(void);
static void test_orient_neg(void);
static void test_orient_robust(void);
static void test_get_trg_area(void);
static void test_trg_get_centroid(void);
static void test_get_circumradius(void);
static void test_get_circumcenter(void);
static void test_get_circumcenter_from_sgm(void);
static void test_get_trg_min_angle(void);
static void test_get_cr2se_ratio(void);
static void test_get_trg_quality(void);
static void test_get_min_trg_edge(void);
static void test_get_max_trg_edge(void);
static void test_get_closest_vtx(void);
static void test_get_closest_vtx_not_ignored(void);
static void test_get_closest_pnt_to_sgm_case_middle(void);
static void test_get_closest_pnt_to_sgm_case_vtx1(void);
static void test_get_closest_pnt_to_sgm_case_vtx2(void);
static void test_are_sgm_intersected_true(void);
static void test_are_sgm_intersected_parallel(void);
static void test_are_sgm_intersected_strictly_false(void);
static void test_are_sgm_intersected_a1(void);
static void test_are_sgm_intersected_a2(void);
static void test_are_sgm_intersected_b1(void);
static void test_are_sgm_intersected_b2(void);
static void test_sgm_intersects_trg_1edge(void);
static void test_sgm_intersects_trg_2edge(void);
static void test_sgm_intersects_trg_false(void);
static void test_sgm_intersects_circle_1x(void);
static void test_sgm_intersects_circle_2x(void);
static void test_sgm_intersects_circle_false(void);
static void test_pnt_lies_on_sgm_true(void);
static void test_pnt_lies_on_sgm_false(void);
static void test_pnt_lies_in_trg_true(void);
static void test_pnt_lies_in_trg_true_on_limit(void);
static void test_pnt_lies_in_trg_false(void);
static void test_pnt_lies_strictly_in_diametral_circle_true(void);
static void test_pnt_lies_strictly_in_diametral_circle_false_on_limit(void);
static void test_pnt_lies_strictly_in_diametral_circle_false(void);
static void test_pnt_lies_strictly_in_circumcircle_true(void);
static void test_pnt_lies_strictly_in_circumcircle_false_on_limit(void);
static void test_pnt_lies_strictly_in_circumcircle_false(void);
static void test_pnt_lies_in_box_true(void);
static void test_pnt_lies_in_box_true_on_limit(void);
static void test_pnt_lies_in_box_false(void);

void cunit_nb_geometric_bot_utils2D(void)
{
	CU_pSuite suite = CU_add_suite("nb/geometric_bot/utils2D.c",
				       suite_init, suite_clean);
	CU_add_test(suite, "get_x_from_darray()",
		    test_get_x_from_darray);
	CU_add_test(suite, "get_y_from_darray()",
		    test_get_y_from_darray);
	CU_add_test(suite, "get_normal()", test_get_normal);
	CU_add_test(suite, "get_dist()", test_get_dist);
	CU_add_test(suite, "get_dist2()", test_get_dist2);
	CU_add_test(suite, "get_delaunay_dist()",
		    test_get_delaunay_dist);
	CU_add_test(suite, "get_delaunay_dist_neg()",
		    test_get_delaunay_dist_neg);
	CU_add_test(suite, "get_enveloping_box()",
		    test_get_enveloping_box);
	CU_add_test(suite, "get_enveloping_box_from_subset()",
		    test_get_enveloping_box_from_subset);
	CU_add_test(suite, "orient()", test_orient);
	CU_add_test(suite, "orient_neg()", test_orient_neg);
	CU_add_test(suite, "orient_robust()", test_orient_robust);
	CU_add_test(suite, "get_trg_area()", test_get_trg_area);
	CU_add_test(suite, "trg_get_centroid()",
		    test_trg_get_centroid);
	CU_add_test(suite, "get_circumradius()",
		    test_get_circumradius);
	CU_add_test(suite, "get_circumcenter()",
		    test_get_circumcenter);
	CU_add_test(suite, "get_circumcenter_from_sgm()",
		    test_get_circumcenter_from_sgm);
	CU_add_test(suite, "get_trg_min_angle()",
		    test_get_trg_min_angle);
	CU_add_test(suite, "get_cr2se_ratio()", test_get_cr2se_ratio);
	CU_add_test(suite, "get_trg_quality()", test_get_trg_quality);
	CU_add_test(suite, "get_min_trg_edge()",
		    test_get_min_trg_edge);
	CU_add_test(suite, "get_max_trg_edge()",
		    test_get_max_trg_edge);
	CU_add_test(suite, "get_closest_vtx()", test_get_closest_vtx);
	CU_add_test(suite, "get_closest_vtx_not_ignored()",
		    test_get_closest_vtx_not_ignored);
	CU_add_test(suite,
		    "get_closest_pnt_to_sgm() if closest is in sgm",
		    test_get_closest_pnt_to_sgm_case_middle);
	CU_add_test(suite,
		    "get_closest_pnt_to_sgm() if closest is vtx 1",
		    test_get_closest_pnt_to_sgm_case_vtx1);
	CU_add_test(suite,
		    "get_closest_pnt_to_sgm() if closest is vtx 2",
		    test_get_closest_pnt_to_sgm_case_vtx2);
	CU_add_test(suite,
		    "are_sgm_intersected() with true intersection",
		    test_are_sgm_intersected_true);
	CU_add_test(suite,
		    "are_sgm_intersected() if sgm are parallel",
		    test_are_sgm_intersected_parallel);
	CU_add_test(suite,
		    "are_sgm_intersected() if sgm are not intersected",
		    test_are_sgm_intersected_strictly_false);
	CU_add_test(suite,
		    "are_sgm_intersected() if a1 is the intersection",
		    test_are_sgm_intersected_a1);
	CU_add_test(suite,
		    "are_sgm_intersected() if a2 is the intersection",
		    test_are_sgm_intersected_a2);
	CU_add_test(suite,
		    "are_sgm_intersected() if b1 is the intersection",
		    test_are_sgm_intersected_b1);
	CU_add_test(suite,
		    "are_sgm_intersected() if b2 is the intersection",
		    test_are_sgm_intersected_b2);
	CU_add_test(suite,
		    "sgm_intersects_trg() with one edge intersection",
		    test_sgm_intersects_trg_1edge);
	CU_add_test(suite,
		    "sgm_intersects_trg() with two edges intersection",
		    test_sgm_intersects_trg_2edge);
	CU_add_test(suite,
		    "sgm_intersects_trg() with no intersections",
		    test_sgm_intersects_trg_false);
	CU_add_test(suite,
		    "sgm_intersects_circle() with one intersection",
		    test_sgm_intersects_circle_1x);
	CU_add_test(suite,
		    "sgm_intersects_circle() with two intersection",
		    test_sgm_intersects_circle_2x);
	CU_add_test(suite,
		    "sgm_intersects_circle() with no intersections",
		    test_sgm_intersects_circle_false);
	CU_add_test(suite,
		    "pnt_lies_on_sgm() if lies inside",
		    test_pnt_lies_on_sgm_true);
	CU_add_test(suite,
		    "pnt_lies_on_sgm() if lies outside",
		    test_pnt_lies_on_sgm_false);
	CU_add_test(suite,
		    "pnt_lies_in_trg() if lies inside",
		    test_pnt_lies_in_trg_true);
	CU_add_test(suite,
		    "pnt_lies_in_trg() if lies in the edge",
		    test_pnt_lies_in_trg_true_on_limit);
	CU_add_test(suite,
		    "pnt_lies_in_trg() if lies outside",
		    test_pnt_lies_in_trg_false);
	CU_add_test(suite,
		    "pnt_lies_strict_in_diam_circ() if lies inside",
		    test_pnt_lies_strictly_in_diametral_circle_true);
	CU_add_test(suite,
		    "pnt_lies_strict_in_diam_circ() if lies on limits",
		    test_pnt_lies_strictly_in_diametral_circle_false_on_limit);
	CU_add_test(suite,
		    "pnt_lies_strict_in_diam_circ() if lies outside",
		    test_pnt_lies_strictly_in_diametral_circle_false);
	CU_add_test(suite,
		    "pnt_lies_strict_in_circumc() if lies inside",
		    test_pnt_lies_strictly_in_circumcircle_true);
	CU_add_test(suite,
		    "pnt_lies_strict_in_circumc() if lies on limits",
		    test_pnt_lies_strictly_in_circumcircle_false_on_limit);
	CU_add_test(suite,
		    "pnt_lies_strict_in_circumc() if lies outside",
		    test_pnt_lies_strictly_in_circumcircle_false);
	CU_add_test(suite,
		    "pnt_lies_in_box() if lies inside",
		    test_pnt_lies_in_box_true);
	CU_add_test(suite,
		    "pnt_lies_in_box() if lies on limits",
		    test_pnt_lies_in_box_true_on_limit);
	CU_add_test(suite,
		    "pnt_lies_in_box() if lies outside",
		    test_pnt_lies_in_box_false);
}

static int suite_init(void)
{
	return 0;
}

static int suite_clean(void)
{
	return 0;
}

static void test_get_x_from_darray(void)
{
	int id = 2;
	double vertices[10] = {0, 5, 2, 1, 8, 9, 6, 3, 7, 4};
	void *vtx = &(vertices[id * 2]);
	double x = vcn_utils2D_get_x_from_darray(vtx);
	CU_ASSERT(vertices[id * 2] == x);
}

static void test_get_y_from_darray(void)
{
	int id = 3;
	double vertices[10] = {0, 5, 2, 1, 8, 9, 6, 3, 7, 4};
	void *vtx = &(vertices[id * 2]);
	double y = vcn_utils2D_get_y_from_darray(vtx);
	CU_ASSERT(vertices[id*2+1] == y);
}

static void test_get_normal(void)
{
	double x1[2] = {0, 0};
	double x2[2] = {1, 1};
	double normal[2];
	vcn_utils2D_get_normal(x1, x2, normal);
	double dist = sqrt(POW2(normal[0]) + POW2(normal[1]));
	CU_ASSERT(fabs(1.0 - dist) < TOLERANCE);
	CU_ASSERT(fabs(normal[0] - INV_SQRT2) < TOLERANCE);
	CU_ASSERT(fabs(normal[1] - INV_SQRT2) < TOLERANCE);
}

static void test_get_dist(void)
{
	double p1[2] = {0, 0};
	double p2[2] = {1, 1};
	double dist = vcn_utils2D_get_dist(p1, p2);
	CU_ASSERT(fabs(dist - SQRT2) < TOLERANCE);
}

static void test_get_dist2(void)
{
	double p1[2] = {0, 0};
	double p2[2] = {1, 1};
	double dist = vcn_utils2D_get_dist2(p1, p2);
	CU_ASSERT(fabs(dist - 2) < TOLERANCE);
}

static void test_get_delaunay_dist(void)
{
	double p1[2] = {0, 0};
	double p2[2] = {1, 0};
	double p3[2] = {1, 1};
	double dist = vcn_utils2D_get_delaunay_dist(p1, p2, p3);
	CU_ASSERT(fabs(dist - 0.5) < TOLERANCE);
}

static void test_get_delaunay_dist_neg(void)
{
	double p1[2] = {0, 0};
	double p2[2] = {2, 0};
	double p3[2] = {1, 0.5};
	double dist = vcn_utils2D_get_delaunay_dist(p1, p2, p3);
	CU_ASSERT(dist < 0.0);
}

static void test_get_enveloping_box(void)
{
	double vertices[10] = {0, 0, 0, 1, 1, 1, 1, 0, 0.5, 0.5};
	int16_t vtx_size = 2 * sizeof(double);
	double box[4];
	vcn_utils2D_get_enveloping_box(5, vertices, vtx_size,
				       vcn_utils2D_get_x_from_darray,
				       vcn_utils2D_get_y_from_darray,
				       box);
	CU_ASSERT(fabs(box[0] - 0) < TOLERANCE);
	CU_ASSERT(fabs(box[1] - 0) < TOLERANCE);
	CU_ASSERT(fabs(box[2] - 1) < TOLERANCE);
	CU_ASSERT(fabs(box[3] - 1) < TOLERANCE);
}

static void test_get_enveloping_box_from_subset(void)
{
	double vertices[10] = {0, 0, 0, 1, 1, 1, 1, 0, 0.5, 0.5};
	int16_t vtx_size = 2 * sizeof(double);
	uint32_t subset[2] = {0, 4};
	double box[4];
	vcn_utils2D_get_enveloping_box_from_subset
		(2, subset,
		 vertices, vtx_size,
		 vcn_utils2D_get_x_from_darray,
		 vcn_utils2D_get_y_from_darray,
		 box);
	CU_ASSERT(fabs(box[0] - 0) < TOLERANCE);
	CU_ASSERT(fabs(box[1] - 0) < TOLERANCE);
	CU_ASSERT(fabs(box[2] - 0.5) < TOLERANCE);
	CU_ASSERT(fabs(box[3] - 0.5) < TOLERANCE);
}

static void test_orient(void)
{
	double t1[2] = {0, 0};
	double t2[2] = {1, 0};
	double t3[2] = {1, 1};
	double orient = vcn_utils2D_orient(t1, t2, t3);
	CU_ASSERT(fabs(orient - 1.0) < TOLERANCE);
}

static void test_orient_neg(void)
{
	double t1[2] = {0, 0};
	double t2[2] = {1, 0};
	double t3[2] = {1, 1};
	double orient = vcn_utils2D_orient(t2, t1, t3);
	CU_ASSERT(fabs(orient + 1.0) < TOLERANCE);
}

static void test_orient_robust(void)
{
	double t1[2] = {1e12, 2e12};
	double t2[2] = {1.1e24, 2.2e24};
	double t3[2] = {1000000.1000001e-12, 2000000.2e-12};
	double orient = vcn_utils2D_orient(t1, t2, t3);
	CU_ASSERT(fabs(orient + 219889.753107016469585) < TOLERANCE);
}

static void test_get_trg_area(void)
{
	double t1[2] = {0, 0};
	double t2[2] = {1, 0};
	double t3[2] = {1, 1};
	double area = vcn_utils2D_get_trg_area(t1, t2, t3);
	CU_ASSERT(fabs(area - 0.5) < TOLERANCE);
}

static void test_trg_get_centroid(void)
{
	double t1[2] = {0, 0};
	double t2[2] = {1, 0};
	double t3[2] = {1, 1};
	double centroid[2];
	vcn_utils2D_trg_get_centroid(t1, t2, t3, centroid);
	CU_ASSERT(fabs(centroid[0] - 2.0/3.0) < TOLERANCE);
	CU_ASSERT(fabs(centroid[1] - 1.0/3.0) < TOLERANCE);
}

static void test_get_circumradius(void)
{
	double t1[2] = {0, 0};
	double t2[2] = {1, 0};
	double t3[2] = {1, 1};
	double circumradius = vcn_utils2D_get_circumradius(t1, t2, t3);
	CU_ASSERT(fabs(circumradius - SQRT2/2.0) < TOLERANCE);
}

static void test_get_circumcenter(void)
{
	double t1[2] = {0, 0};
	double t2[2] = {1, 0};
	double t3[2] = {1, 1};
	double circumcenter[2];
	vcn_utils2D_get_circumcenter(t1, t2, t3, circumcenter);
	CU_ASSERT(fabs(circumcenter[0] - 0.5) < TOLERANCE);
	CU_ASSERT(fabs(circumcenter[1] - 0.5) < TOLERANCE);
}

static void test_get_circumcenter_from_sgm(void)
{
	double s1[2] = {0, 0};
	double s2[2] = {1, 0};
	double circumcenter[2];
	vcn_utils2D_get_circumcenter_from_sgm(s1, s2, SQRT2/2.0,
					      circumcenter);
	CU_ASSERT(fabs(circumcenter[0] - 0.5) < TOLERANCE);
	CU_ASSERT(fabs(circumcenter[1] - 0.5) < TOLERANCE);
}

static void test_get_trg_min_angle(void)
{
	double t1[2] = {0, 0};
	double t2[2] = {1, 0};
	double t3[2] = {1, 1};
	double min_angle = vcn_utils2D_get_trg_min_angle(t1, t2, t3);
	CU_ASSERT(fabs(min_angle - PI/4.0) < TOLERANCE);
}

static void test_get_cr2se_ratio(void)
{
	double t1[2] = {0, 0};
	double t2[2] = {1, 0};
	double t3[2] = {1, 1};
	double cr2se = vcn_utils2D_get_cr2se_ratio(t1, t2, t3);
	CU_ASSERT(fabs(cr2se - SQRT2/2) < TOLERANCE);
}

static void test_get_trg_quality(void)
{
	double t1[2] = {0, 0};
	double t2[2] = {1, 0};
	double t3[2] = {1, 1};
	double quality = vcn_utils2D_get_trg_quality(t1, t2, t3);
	double cr2se = vcn_utils2D_get_cr2se_ratio(t1, t2, t3);
	CU_ASSERT(fabs(quality - INV_SQRT3/cr2se) < TOLERANCE);
}

static void test_get_min_trg_edge(void)
{
	double t1[2] = {0, 0};
	double t2[2] = {1, 0};
	double t3[2] = {1, 1};
	double min = vcn_utils2D_get_min_trg_edge(t1, t2, t3);
	CU_ASSERT(fabs(min - 1.0) < TOLERANCE);
}

static void test_get_max_trg_edge(void)
{
	double t1[2] = {0, 0};
	double t2[2] = {1, 0};
	double t3[2] = {1, 1};
	double max = vcn_utils2D_get_max_trg_edge(t1, t2, t3);
	CU_ASSERT(fabs(max - SQRT2) < TOLERANCE);
}

static void test_get_closest_vtx(void)
{
	double vertices[10] = {0, 0, 0, 1, 1, 1, 1, 0, 0.5, 0.5};
	uint32_t id_closest =
		vcn_utils2D_get_closest_vtx(2, 2, 5, vertices);
	CU_ASSERT(2 == id_closest);

}

static void test_get_closest_vtx_not_ignored(void)
{
	double vertices[10] = {0, 0, 0, 1, 1, 1, 1, 0, 0.5, 0.5};
	uint32_t ignore[2] = {0, 2};
	uint32_t id_closest =
		vcn_utils2D_get_closest_vtx_not_ignored(-1, -1,
							5, vertices,
							2, ignore);
	CU_ASSERT(4 == id_closest);

}

static void test_get_closest_pnt_to_sgm_case_middle(void)
{
	double s1[2] = {0, 0};
	double s2[2] = {1, 1};
	double p[2] = {0, 1};
	double cp[2];
	vcn_utils2D_get_closest_pnt_to_sgm(s1, s2, p, cp);
	CU_ASSERT(fabs(cp[0] - 0.5) < TOLERANCE);
	CU_ASSERT(fabs(cp[1] - 0.5) < TOLERANCE);
}

static void test_get_closest_pnt_to_sgm_case_vtx1(void)
{
	double s1[2] = {0, 0};
	double s2[2] = {1, 1};
	double p[2] = {0, -1};
	double cp[2];
	vcn_utils2D_get_closest_pnt_to_sgm(s1, s2, p, cp);
	CU_ASSERT(fabs(cp[0] - 0) < TOLERANCE);
	CU_ASSERT(fabs(cp[1] - 0) < TOLERANCE);
}

static void test_get_closest_pnt_to_sgm_case_vtx2(void)
{
	double s1[2] = {0, 0};
	double s2[2] = {1, 1};
	double p[2] = {1, 2};
	double cp[2];
	vcn_utils2D_get_closest_pnt_to_sgm(s1, s2, p, cp);
	CU_ASSERT(fabs(cp[0] - 1) < TOLERANCE);
	CU_ASSERT(fabs(cp[1] - 1) < TOLERANCE);
}

static void test_are_sgm_intersected_true(void)
{
	double a1[2] = {0, 0};
	double a2[2] = {1, 1};
	double b1[2] = {1, 0};
	double b2[2] = {0, 1};
	double intersection[2];
	nb_intersect_t status =
		vcn_utils2D_get_sgm_intersection(a1, a2, b1, b2,
						 intersection);
	CU_ASSERT(fabs(intersection[0] - 0.5) < TOLERANCE);
	CU_ASSERT(fabs(intersection[1] - 0.5) < TOLERANCE);
	CU_ASSERT(NB_INTERSECTED == status);
}

static void test_are_sgm_intersected_parallel(void)
{
	double a1[2] = {0, 0};
	double a2[2] = {0, 1};
	double b1[2] = {1, 0};
	double b2[2] = {1, 1};
	double intersection[2];
	nb_intersect_t status =
		vcn_utils2D_get_sgm_intersection(a1, a2, b1, b2,
						 intersection);
	CU_ASSERT(NB_PARALLEL == status);
}

static void test_are_sgm_intersected_strictly_false(void)
{
	double a1[2] = {0, 0};
	double a2[2] = {1, 1};
	double b1[2] = {2, 1};
	double b2[2] = {1, 2};
	double intersection[2];
	nb_intersect_t status =
		vcn_utils2D_get_sgm_intersection(a1, a2, b1, b2,
						 intersection);
	CU_ASSERT(NB_NOT_INTERSECTED == status);
}

static void test_are_sgm_intersected_a1(void)
{
	double a1[2] = {0, 0};
	double a2[2] = {1, 1};
	double b1[2] = {-1, 1};
	double b2[2] = {1, -1};
	double intersection[2];
	nb_intersect_t status =
		vcn_utils2D_get_sgm_intersection(a1, a2, b1, b2,
						 intersection);
	CU_ASSERT(fabs(intersection[0] - a1[0]) < TOLERANCE);
	CU_ASSERT(fabs(intersection[1] - a1[1]) < TOLERANCE);
	CU_ASSERT(NB_INTERSECT_ON_A1 == status);
}

static void test_are_sgm_intersected_a2(void)
{
	double a1[2] = {1, 1};
	double a2[2] = {0, 0};
	double b1[2] = {-1, 1};
	double b2[2] = {1, -1};
	double intersection[2];
	nb_intersect_t status =
		vcn_utils2D_get_sgm_intersection(a1, a2, b1, b2,
						 intersection);
	CU_ASSERT(fabs(intersection[0] - a2[0]) < TOLERANCE);
	CU_ASSERT(fabs(intersection[1] - a2[1]) < TOLERANCE);
	CU_ASSERT(NB_INTERSECT_ON_A2 == status);
}

static void test_are_sgm_intersected_b1(void)
{
	double a1[2] = {-1, 1};
	double a2[2] = {1, -1};
	double b1[2] = {0, 0};
	double b2[2] = {1, 1};
	double intersection[2];
	nb_intersect_t status =
		vcn_utils2D_get_sgm_intersection(a1, a2, b1, b2,
						 intersection);
	CU_ASSERT(fabs(intersection[0] - b1[0]) < TOLERANCE);
	CU_ASSERT(fabs(intersection[1] - b1[1]) < TOLERANCE);
	CU_ASSERT(NB_INTERSECT_ON_B1 == status);
}

static void test_are_sgm_intersected_b2(void)
{
	double a1[2] = {-1, 1};
	double a2[2] = {1, -1};
	double b1[2] = {1, 1};
	double b2[2] = {0, 0};
	double intersection[2];
	nb_intersect_t status =
		vcn_utils2D_get_sgm_intersection(a1, a2, b1, b2,
						 intersection);
	CU_ASSERT(fabs(intersection[0] - b2[0]) < TOLERANCE);
	CU_ASSERT(fabs(intersection[1] - b2[1]) < TOLERANCE);
	CU_ASSERT(NB_INTERSECT_ON_B2 == status);
}

static void test_sgm_intersects_trg_1edge(void)
{
	double t1[2] = {0, 0};
	double t2[2] = {1, 0};
	double t3[2] = {1, 1};
	double s1[2] = {1, 0};
	double s2[2] = {0.2, 0.8};
	bool intersects =
		vcn_utils2D_sgm_intersects_trg(t1, t2, t3, s1, s2);
	CU_ASSERT(intersects);
}

static void test_sgm_intersects_trg_2edge(void)
{
	double t1[2] = {0, 0};
	double t2[2] = {1, 0};
	double t3[2] = {1, 1};
	double s1[2] = {1, 0};
	double s2[2] = {0.2, 1.2};
	bool intersects =
		vcn_utils2D_sgm_intersects_trg(t1, t2, t3, s1, s2);
	CU_ASSERT(intersects);
}

static void test_sgm_intersects_trg_false(void)
{
	double t1[2] = {0, 0};
	double t2[2] = {1, 0};
	double t3[2] = {1, 1};
	double s1[2] = {1, 0};
	double s2[2] = {0.2, -0.8};
	bool not_intersects = 
		!vcn_utils2D_sgm_intersects_trg(t1, t2, t3, s1, s2);
	CU_ASSERT(not_intersects);
}

static void test_sgm_intersects_circle_1x(void)
{
	double s1[2] = {0, 0};
	double s2[2] = {1, 1};
	double cc[2] = {1, 1};
	bool intersects = 
		vcn_utils2D_sgm_intersects_circle(cc, 0.5, s1, s2);
	CU_ASSERT(intersects);
}

static void test_sgm_intersects_circle_2x(void)
{
	double s1[2] = {0, 0};
	double s2[2] = {2, 2};
	double cc[2] = {1, 1};
	bool intersects = vcn_utils2D_sgm_intersects_circle(cc, 0.5, s1, s2);
	CU_ASSERT(intersects);
}

static void test_sgm_intersects_circle_false(void)
{
	double s1[2] = {0, 0};
	double s2[2] = {-1, -1};
	double cc[2] = {1, 1};
	bool not_intersects = 
		!vcn_utils2D_sgm_intersects_circle(cc, 0.5, s1, s2);
	CU_ASSERT(not_intersects);
}

static void test_pnt_lies_on_sgm_true(void)
{
	double s1[2] = {0, 0};
	double s2[2] = {1, 1};
	double p[2] = {0.5, 0.5};
	bool lies_inside = vcn_utils2D_pnt_lies_on_sgm(s1, s2, p);
	CU_ASSERT(lies_inside);
}

static void test_pnt_lies_on_sgm_false(void)
{
	double s1[2] = {0, 0};
	double s2[2] = {1, 1};
	double p[2] = {1.1, 1.1};
	bool lies_outside = !vcn_utils2D_pnt_lies_on_sgm(s1, s2, p);
	CU_ASSERT(lies_outside);
}

static void test_pnt_lies_in_trg_true(void)
{
	double t1[2] = {0, 0};
	double t2[2] = {1, 0};
	double t3[2] = {1, 1};
	double p[2] = {0.6, 0.5};
	bool lies_inside = vcn_utils2D_pnt_lies_in_trg(t1, t2, t3, p);
	CU_ASSERT(lies_inside);
}

static void test_pnt_lies_in_trg_true_on_limit(void)
{
	double t1[2] = {0, 0};
	double t2[2] = {1, 0};
	double t3[2] = {1, 1};
	double p[2] = {0.5, 0.5};
	bool lies_inside = vcn_utils2D_pnt_lies_in_trg(t1, t2, t3, p);
	CU_ASSERT(lies_inside);
}

static void test_pnt_lies_in_trg_false(void)
{
	double t1[2] = {0, 0};
	double t2[2] = {1, 0};
	double t3[2] = {1, 1};
	double p[2] = {0.4, 0.5};
	bool lies_outside =
		!vcn_utils2D_pnt_lies_in_trg(t1, t2, t3, p);
	CU_ASSERT(lies_outside);
}

static void test_pnt_lies_strictly_in_diametral_circle_true(void)
{
	double s1[2] = {0, 0};
	double s2[2] = {1, 0};
	double p[2] = {0.5, 0};
	bool lies_inside =
		vcn_utils2D_pnt_lies_strictly_in_diametral_circle(s1, s2, p);
	CU_ASSERT(lies_inside);
}

static void test_pnt_lies_strictly_in_diametral_circle_false_on_limit(void)
{
	double s1[2] = {0, 0};
	double s2[2] = {1, 0};
	double p[2] = {0.5, 0.5};
	bool lies_outside =
		!vcn_utils2D_pnt_lies_strictly_in_diametral_circle(s1, s2, p);
	CU_ASSERT(lies_outside);
}

static void test_pnt_lies_strictly_in_diametral_circle_false(void)
{
	double s1[2] = {0, 0};
	double s2[2] = {1, 0};
	double p[2] = {1, 0.5};
	bool lies_outside =
		!vcn_utils2D_pnt_lies_strictly_in_diametral_circle(s1, s2, p);
	CU_ASSERT(lies_outside);
}

static void test_pnt_lies_strictly_in_circumcircle_true(void)
{
	double t1[2] = {0, 0};
	double t2[2] = {1, 0};
	double t3[2] = {1, 1};
	double p[2] = {0.5, 0.5};
	bool lies_inside =
		vcn_utils2D_pnt_lies_strictly_in_circumcircle(t1, t2, t3, p);
	CU_ASSERT(lies_inside);
}

static void test_pnt_lies_strictly_in_circumcircle_false_on_limit(void)
{
	double t1[2] = {0, 0};
	double t2[2] = {1, 0};
	double t3[2] = {1, 1};
	double p[2] = {1, 0};
	bool lies_outside =
		!vcn_utils2D_pnt_lies_strictly_in_circumcircle(t1, t2, t3, p);
	CU_ASSERT(lies_outside);
}

static void test_pnt_lies_strictly_in_circumcircle_false(void)
{
	double t1[2] = {0, 0};
	double t2[2] = {1, 0};
	double t3[2] = {1, 1};
	double p[2] = {2, 0};
	bool lies_outside =
		!vcn_utils2D_pnt_lies_strictly_in_circumcircle(t1, t2, t3, p);
	CU_ASSERT(lies_outside);
}

static void test_pnt_lies_in_box_true(void)
{
	double box[4] = {0, 0, 1, 1};
	double p[2] = {0.5, 0.5};
	bool lies_inside = vcn_utils2D_pnt_lies_in_box(box, p);
	CU_ASSERT(lies_inside);
}

static void test_pnt_lies_in_box_true_on_limit(void)
{
	double box[4] = {0, 0, 1, 1};
	double p[2] = {0.5, 1};
	bool lies_inside = vcn_utils2D_pnt_lies_in_box(box, p);
	CU_ASSERT(lies_inside);
}

static void test_pnt_lies_in_box_false(void)
{
	double box[4] = {0, 0, 1, 1};
	double p[2] = {0.5, 1.1};
	bool lies_outside = !vcn_utils2D_pnt_lies_in_box(box, p);
	CU_ASSERT(lies_outside);
}
