#include <stdbool.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "nb/container_bot.h"
#include "nb/geometric_bot/utils2D.h"
#include "nb/geometric_bot/knn/bins2D.h"
#include "nb/geometric_bot/knn/bins2D_iterator.h"

#include "test_library.h"
#include "test_add.h"

#define POW2(a) ((a)*(a))

static bool check_create(void);
static bool check_destroy(void);
static bool check_clear(void);
static bool check_enable_point_destroyer(void);
static bool check_disable_point_destroyer(void);
static bool check_set_attribute_destroyer(void);
static bool check_insert(void);
static bool check_delete(void);
static bool check_delete_first(void);
static bool check_get_knn(void);
static bool check_set_filter(void);
static bool check_set_filter_data(void);
static bool check_get_candidate_points_to_min_delaunay(void);
static bool check_get_points_inside_circle(void);
static bool check_are_points_inside_circle(void);
static bool check_get_N_bins(void);
static bool check_get_min_points_x_bin(void);
static bool check_get_length(void);
static bool check_is_empty(void);
static bool check_is_not_empty(void);
static bool check_get_size_of_bins(void);

static vcn_bins2D_t* get_bins(int N);
static vcn_bins2D_t* get_bins_and_array(int N, vcn_point2D_t *vertices[]);
static vcn_bins2D_t* get_bins_with_attribute(int N);
static bool distances_are_ok(int N, double dist[]);
static bool filter(const vcn_point2D_t *const p_ref,
		   const vcn_point2D_t *const p,
		   const void *const data);
static bool all_pass_filter(int N, vcn_point2D_t *p_ref,
			    vcn_point2D_t *knn[]);
static bool filter_with_data(const vcn_point2D_t *const p_ref,
			     const vcn_point2D_t *const p,
			     const void *const data);
static bool contains_min_dd(const vcn_bins2D_t *const bins,
			       const nb_container_t *const cnt,
			       const vcn_point2D_t *const p1,
			       const vcn_point2D_t *const p2);
static bool all_in_half_space(const nb_container_t *const cnt,
			      const vcn_point2D_t *const p1,
			      const vcn_point2D_t *const p2);
static bool dont_include_edge_points(const nb_container_t *const cnt,
				     const vcn_point2D_t *const p1,
				     const vcn_point2D_t *const p2);
static double get_min_dd_in_container(const nb_container_t *const cnt,
				      const vcn_point2D_t *const p1,
				      const vcn_point2D_t *const p2);
static double get_min_dd_in_bins(const vcn_bins2D_t *const bins,
				 const vcn_point2D_t *const p1,
				 const vcn_point2D_t *const p2);
static bool all_inside_circle(const nb_container_t *const cnt,
			      const double center[2],
			      double radius);

inline int vcn_test_get_driver_id(void)
{
	return NB_DRIVER_UNIT_TEST;
}

void vcn_test_load_tests(void *tests_ptr)
{
	vcn_test_add(tests_ptr, check_create,
		     "Check create()");
	vcn_test_add(tests_ptr, check_destroy,
		     "Check destroy()");
	vcn_test_add(tests_ptr, check_clear,
		     "Check clear()");
	vcn_test_add(tests_ptr, check_enable_point_destroyer,
		     "Check enable_point_destroyer()");
	vcn_test_add(tests_ptr, check_disable_point_destroyer,
		     "Check disable_point_destroyer()");
	vcn_test_add(tests_ptr, check_set_attribute_destroyer,
		     "Check set_attribute_destroyer()");
	vcn_test_add(tests_ptr, check_insert,
		     "Check insert()");
	vcn_test_add(tests_ptr, check_delete,
		     "Check delete()");
	vcn_test_add(tests_ptr, check_delete_first,
		     "Check delete_first()");
	vcn_test_add(tests_ptr, check_get_knn,
		     "Check get_knn()");
	vcn_test_add(tests_ptr, check_set_filter,
		     "Check set_filter()");
	vcn_test_add(tests_ptr, check_set_filter_data,
		     "Check set_filter_data()");
	vcn_test_add(tests_ptr, check_get_candidate_points_to_min_delaunay,
		     "Check get_candidate_points_to_min_delaunay()");
	vcn_test_add(tests_ptr, check_get_points_inside_circle,
		     "Check get_points_inside_circle()");
	vcn_test_add(tests_ptr, check_are_points_inside_circle,
		     "Check are_points_inside_circle()");
	vcn_test_add(tests_ptr, check_get_N_bins,
		     "Check get_N_bins()");
	vcn_test_add(tests_ptr, check_get_min_points_x_bin,
		     "Check get_min_points_x_bin()");
	vcn_test_add(tests_ptr, check_get_length,
		     "Check get_length()");
	vcn_test_add(tests_ptr, check_is_empty,
		     "Check is_empty()");
	vcn_test_add(tests_ptr, check_is_not_empty,
		     "Check is_not_empty()");
	vcn_test_add(tests_ptr, check_get_size_of_bins,
		     "Check get_size_of_bins()");
}

static bool check_create(void)
{
	vcn_bins2D_t *bins = vcn_bins2D_create(1.0);
	bool is_ok = (NULL != bins);
	if (is_ok) {
		bool size_is_ok = 
			(fabs(1.0 - vcn_bins2D_get_size_of_bins(bins)) < 1e-9);
		bool length_is_ok = (0 == vcn_bins2D_get_length(bins));
		is_ok = is_ok && size_is_ok && length_is_ok;
		vcn_bins2D_destroy(bins);
	}
	return is_ok;
}

static bool check_destroy(void)
{
	vcn_bins2D_t *bins = vcn_bins2D_create(1.0);
	vcn_bins2D_destroy(bins);
	return true;
}

static bool check_clear(void)
{
	int N = 100;
	vcn_bins2D_t *bins = get_bins(N);
	vcn_bins2D_clear(bins);
	bool is_ok = (0 == vcn_bins2D_get_length(bins));
	vcn_bins2D_destroy(bins);
	return is_ok;
}

static bool check_enable_point_destroyer(void)
{
	int N = 100;
	vcn_bins2D_t *bins = get_bins(N);
	vcn_bins2D_disable_point_destroyer(bins);
	vcn_bins2D_enable_point_destroyer(bins);
	vcn_bins2D_destroy(bins);
	return true;
}

static bool check_disable_point_destroyer(void)
{
	int N = 100;
	vcn_point2D_t *vertices[100];
	vcn_bins2D_t *bins = get_bins_and_array(N, vertices);
	vcn_bins2D_disable_point_destroyer(bins);
	vcn_bins2D_destroy(bins);
	for (int i = 0; i < N; i++)
		vcn_point2D_destroy(vertices[i]);
	return true;
}

static bool check_set_attribute_destroyer(void)
{
	int N = 100;
	vcn_bins2D_t *bins = get_bins_with_attribute(N);
	vcn_bins2D_set_attribute_destroyer(bins, free);
	vcn_bins2D_destroy(bins);
	return true;
}

static bool check_insert(void)
{
	int N = 100;
	vcn_bins2D_t *bins = get_bins(N);
	bool is_ok = (N == vcn_bins2D_get_length(bins));
	vcn_bins2D_destroy(bins);
	return is_ok;
}

static bool check_delete(void)
{
	int N = 100;
	vcn_bins2D_t *bins = get_bins(N);
	vcn_point2D_t point;
	point.x[0] = 0.0;
	point.x[1] = 0.0;
	vcn_point2D_t *deleted_point = vcn_bins2D_delete(bins, &point);
	bool is_ok = false;
	if (NULL != deleted_point) {
		is_ok = (fabs(deleted_point->x[0] - point.x[0]) < 1e-9) &&
			(fabs(deleted_point->x[1] - point.x[1]) < 1e-9);
		vcn_point2D_destroy(deleted_point);
	}
	vcn_bins2D_destroy(bins);
	return is_ok;
}

static bool check_delete_first(void)
{
	int N = 100;
	vcn_bins2D_t *bins = get_bins(N);
	vcn_point2D_t *deleted_point = vcn_bins2D_delete_first(bins);
	bool is_ok = false;
	if (NULL != deleted_point) {
		vcn_point2D_destroy(deleted_point);
		is_ok = true;
	}
	vcn_bins2D_destroy(bins);
	return is_ok;
}

static bool check_get_knn(void)
{
	int N = 100;
	vcn_bins2D_t *bins = get_bins(N);
	vcn_point2D_t p;
	p.x[0] = 0.5;
	p.x[1] = 0.5;
	vcn_point2D_t* knn[10];
	double knn_dist[10];
	uint32_t k = vcn_bins2D_get_knn(bins, &p, 10, knn, knn_dist);
	bool is_ok = (10 == k) && distances_are_ok(10, knn_dist);
	vcn_bins2D_destroy(bins);
	return is_ok;
}

static bool check_set_filter(void)
{
	int N = 100;
	vcn_bins2D_t *bins = get_bins(N);
	vcn_bins2D_set_filter(bins, filter);
	vcn_point2D_t p;
	p.x[0] = 0.5;
	p.x[1] = 0.5;
	vcn_point2D_t* knn[10];
	double knn_dist[10];
	uint32_t k = vcn_bins2D_get_knn(bins, &p, 10, knn, knn_dist);
	bool is_ok = (10 == k) && distances_are_ok(10, knn_dist) &&
		all_pass_filter(10, &p, knn);
	vcn_bins2D_destroy(bins);
	return is_ok;
}

static bool check_set_filter_data(void)
{
	int N = 100;
	vcn_bins2D_t *bins = get_bins(N);
	vcn_bins2D_set_filter(bins, filter_with_data);
	int fdata = 5;
	vcn_bins2D_set_filter_data(bins, &fdata);
	vcn_point2D_t p;
	p.x[0] = 0.5;
	p.x[1] = 0.5;
	vcn_point2D_t* knn[10];
	double knn_dist[10];
	uint32_t k = vcn_bins2D_get_knn(bins, &p, 10, knn, knn_dist);
	bool is_ok = (10 == k) && distances_are_ok(10, knn_dist) &&
		all_pass_filter(10, &p, knn);
	vcn_bins2D_destroy(bins);
	return is_ok;
}

static bool check_get_candidate_points_to_min_delaunay(void)
{
	int N = 100;
	vcn_bins2D_t *bins = get_bins(N);
	vcn_point2D_t p1, p2;
	p1.x[0] = -1;
	p1.x[1] = 1;
	p2.x[0] = 1;
	p2.x[1] = -1;
	nb_container_t *cnt = 
		vcn_bins2D_get_candidate_points_to_min_delaunay(bins, &p1, &p2);
	bool is_ok = nb_container_is_not_empty(cnt);
	if (is_ok)
		is_ok = is_ok && dont_include_edge_points(cnt, &p1, &p2);
	if (is_ok)
		is_ok = is_ok && all_in_half_space(cnt, &p1, &p2);
	if (is_ok)
		is_ok = is_ok && contains_min_dd(bins, cnt, &p1, &p2);
	nb_container_destroy(cnt);
	vcn_bins2D_destroy(bins);
	return is_ok;	
}

static bool check_get_points_inside_circle(void)
{
	int N = 100;
	vcn_bins2D_t *bins = get_bins(N);
	double center[2] = {0, 0};
	double radius = 10.0;
	nb_container_t *cnt = 
		vcn_bins2D_get_points_inside_circle(bins, center, radius);
	bool is_ok = nb_container_is_not_empty(cnt) &&
		all_inside_circle(cnt, center, radius);
	nb_container_destroy(cnt);
	vcn_bins2D_destroy(bins);
	return is_ok;
}

static bool check_are_points_inside_circle(void)
{
	int N = 100;
	vcn_bins2D_t *bins = get_bins(N);
	double center[2] = {0, 0};
	double radius = 10.0;
	bool is_ok = vcn_bins2D_are_points_inside_circle(bins, center, radius);
	vcn_bins2D_destroy(bins);
	return is_ok;
}

static bool check_get_N_bins(void)
{
	int N = 100;
	vcn_bins2D_t *bins = get_bins(N);
	uint32_t N_bins = vcn_bins2D_get_N_bins(bins);
	bool is_ok = (N_bins > 0 && N_bins < N);
	vcn_bins2D_destroy(bins);
	return is_ok;
}

static bool check_get_min_points_x_bin(void)
{
	int N = 100;
	vcn_bins2D_t *bins = get_bins(N);
	bool is_ok = (vcn_bins2D_get_min_points_x_bin(bins) == 1);
	vcn_bins2D_destroy(bins);
	return is_ok;
}

static bool check_get_length(void)
{
	int N = 100;
	vcn_bins2D_t *bins = get_bins(N);
	bool is_ok = (vcn_bins2D_get_length(bins) == N);
	vcn_bins2D_destroy(bins);
	return is_ok;
}

static bool check_is_empty(void)
{
	int N = 100;
	vcn_bins2D_t *bins = get_bins(N);
	bool is_ok = !vcn_bins2D_is_empty(bins);
	vcn_bins2D_destroy(bins);
	return is_ok;
}

static bool check_is_not_empty(void)
{
	int N = 100;
	vcn_bins2D_t *bins = get_bins(N);
	bool is_ok = vcn_bins2D_is_not_empty(bins);
	vcn_bins2D_destroy(bins);
	return is_ok;
}

static bool check_get_size_of_bins(void)
{
	int N = 100;
	vcn_bins2D_t *bins = get_bins(N);
	bool is_ok = (fabs(vcn_bins2D_get_size_of_bins(bins) - 1.0) < 1e-9);
	vcn_bins2D_destroy(bins);
	return is_ok;
}

static vcn_bins2D_t* get_bins(int N)
{
	vcn_bins2D_t *bins = vcn_bins2D_create(1.0);
	for (int i = 0; i < N; i++) {
		vcn_point2D_t *point = vcn_point2D_create();
		double xsign = (0 == (i/2 + 1) % 2)?(1.0):(-1.0);
		double ysign = (1 == (i/2) % 2)?(1.0):(-1.0);
		point->x[0] = i * 0.5 * xsign;
		point->x[1] = i * 0.5 * ysign;
		vcn_bins2D_insert(bins, point);
	}
	return bins;
}

static vcn_bins2D_t* get_bins_and_array(int N, vcn_point2D_t *vertices[])
{
	vcn_bins2D_t *bins = vcn_bins2D_create(1.0);
	for (int i = 0; i < N; i++) {
		vcn_point2D_t *point = vcn_point2D_create();
		double xsign = (0 == (i/2 + 1) % 2)?(1.0):(-1.0);
		double ysign = (1 == (i/2) % 2)?(1.0):(-1.0);
		point->x[0] = i * 0.5 * xsign;
		point->x[1] = i * 0.5 * ysign;
		vcn_bins2D_insert(bins, point);
		vertices[i] = point;
	}
	return bins;
}

static vcn_bins2D_t* get_bins_with_attribute(int N)
{
	vcn_bins2D_t *bins = vcn_bins2D_create(1.0);
	for (int i = 0; i < N; i++) {
		vcn_point2D_t *point = vcn_point2D_create();
		point->x[0] = i;
		point->x[1] = i;
		point->attr = malloc(1);
		vcn_bins2D_insert(bins, point);
	}
	return bins;
}

static inline bool distances_are_ok(int N, double dist[])
{
	bool is_ok = true;
	for (int i = 1; i < N && is_ok; i++)
		is_ok = is_ok && (dist[i-1] <= dist[i]);
	return is_ok;
}

static inline bool filter(const vcn_point2D_t *const p_ref,
			  const vcn_point2D_t *const p,
			  const void *const data)
{
	return (vcn_utils2D_get_dist(p_ref->x, p->x) > 5.0);
}

static inline bool all_pass_filter(int N, vcn_point2D_t *p_ref,
				   vcn_point2D_t *knn[])
{
	bool is_ok = true;
	for (int i = 0; i < N && is_ok; i++)
		is_ok = is_ok && filter(p_ref, knn[i], NULL);
	return is_ok;
}

static inline bool filter_with_data(const vcn_point2D_t *const p_ref,
				    const vcn_point2D_t *const p,
				    const void *const data)
{
	return (vcn_utils2D_get_dist(p_ref->x, p->x) > *((int*)data));
}

static bool all_in_half_space(const nb_container_t *const cnt,
			      const vcn_point2D_t *const p1,
			      const vcn_point2D_t *const p2)
{
	bool is_ok = true;
	nb_iterator_t *iter = nb_iterator_create();
	nb_iterator_set_container(iter, cnt);
	while (nb_iterator_has_more(iter)) {
		const vcn_point2D_t *p3 = nb_iterator_get_next(iter);
		double area = vcn_utils2D_get_2x_trg_area(p1->x, p2->x, p3->x);
		if (0 > area) {
			is_ok = false;
			break;
		}
	}
	nb_iterator_destroy(iter);
	return is_ok;
}

static bool dont_include_edge_points(const nb_container_t *const cnt,
				     const vcn_point2D_t *const p1,
				     const vcn_point2D_t *const p2)
{
	bool is_ok = true;
	nb_iterator_t *iter = nb_iterator_create();
	nb_iterator_set_container(iter, cnt);
	while (nb_iterator_has_more(iter)) {
		const vcn_point2D_t *p3 = nb_iterator_get_next(iter);
		if (p1 == p3 || p2 == p3) {
			is_ok = false;
			break;
		}
	}
	nb_iterator_destroy(iter);
	return is_ok;
}

static bool contains_min_dd(const vcn_bins2D_t *const bins,
			    const nb_container_t *const cnt,
			    const vcn_point2D_t *const p1,
			    const vcn_point2D_t *const p2)
{
        double min_dd = get_min_dd_in_container(cnt, p1, p2);
	double min_dd_in_bins = get_min_dd_in_bins(bins, p1, p2);
	return min_dd <= min_dd_in_bins;
}

static double get_min_dd_in_container(const nb_container_t *const cnt,
				      const vcn_point2D_t *const p1,
				      const vcn_point2D_t *const p2)
{
	double min_dd = 1e30;
	nb_iterator_t *iter = nb_iterator_create();
	nb_iterator_set_container(iter, cnt);
	while (nb_iterator_has_more(iter)) {
		const vcn_point2D_t *p3 = nb_iterator_get_next(iter);
		double min = vcn_utils2D_get_delaunay_dist(p1->x, p2->x,
							   p3->x);
		if (min < min_dd)
			min_dd = min;
	}
	nb_iterator_destroy(iter);
	return min_dd;
}

static double get_min_dd_in_bins(const vcn_bins2D_t *const bins,
				 const vcn_point2D_t *const p1,
				 const vcn_point2D_t *const p2)
{
	double min_dd = 1e30;
	vcn_bins2D_iter_t *iter = vcn_bins2D_iter_create();
	vcn_bins2D_iter_set_bins(iter, bins);
	while (vcn_bins2D_iter_has_more(iter)) {
		const vcn_point2D_t *p3 = vcn_bins2D_iter_get_next(iter);
		if (p1 != p3 && p2 != p3) {
			double area = vcn_utils2D_get_2x_trg_area(p1->x,
								  p2->x,
								  p3->x);
			if (area > 0.0) {
				double min = 
					vcn_utils2D_get_delaunay_dist(p1->x,
								      p2->x,
								      p3->x);
				if (min < min_dd)
					min_dd = min;
			}
		}
	}
	vcn_bins2D_iter_destroy(iter);
	return min_dd;
}

static bool all_inside_circle(const nb_container_t *const cnt,
			      const double center[2],
			      double radius)
{
	nb_iterator_t *iter = nb_iterator_create();
	nb_iterator_set_container(iter, cnt);
	bool all_inside = true;
	while (nb_iterator_has_more(iter)) {
		const vcn_point2D_t *p = nb_iterator_get_next(iter);
		double dist2 = vcn_utils2D_get_dist2(center, p->x);
		if (POW2(radius) - dist2 < 0) {
			all_inside = false;
			break;
		}		
	}
	nb_iterator_destroy(iter);
	return all_inside;
}
