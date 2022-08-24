#include <stdbool.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "cunit/Basic.h"

#include "nb/memory_bot.h"
#include "nb/container_bot.h"
#include "nb/geometric_bot/utils2D.h"
#include "nb/geometric_bot/knn/bins2D.h"
#include "nb/geometric_bot/knn/bins2D_iterator.h"

#define POW2(a) ((a)*(a))

static int suite_init(void);
static int suite_clean(void);

static void test_create(void);
static void test_destroy(void);
static void test_clear(void);
/* static void test_enable_point_destroyer(void); */
/* static void test_disable_point_destroyer(void); */
/* static void test_set_attribute_destroyer(void); */
static void test_insert(void);
static void test_delete(void);
static void test_delete_first(void);
static void test_get_knn(void);
static void test_set_filter(void);
static void test_set_filter_data(void);
static void test_get_candidate_points_to_min_delaunay(void);
static void test_get_points_inside_circle(void);
static void test_are_points_inside_circle(void);
static void test_get_N_bins(void);
static void test_get_min_points_x_bin(void);
static void test_get_length(void);
static void test_is_empty(void);
static void test_is_not_empty(void);
static void test_get_size_of_bins(void);

static nb_bins2D_t* get_bins(int N);
static bool distances_are_ok(int N, double dist[]);
static bool filter(const nb_point2D_t *const p_ref,
		   const nb_point2D_t *const p,
		   const void *const data);
static bool all_pass_filter(int N, nb_point2D_t *p_ref,
			    nb_point2D_t *knn[]);
static bool filter_with_data(const nb_point2D_t *const p_ref,
			     const nb_point2D_t *const p,
			     const void *const data);
static bool contains_min_dd(const nb_bins2D_t *const bins,
			       const nb_container_t *const cnt,
			       const nb_point2D_t *const p1,
			       const nb_point2D_t *const p2);
static bool all_in_half_space(const nb_container_t *const cnt,
			      const nb_point2D_t *const p1,
			      const nb_point2D_t *const p2);
static bool dont_include_edge_points(const nb_container_t *const cnt,
				     const nb_point2D_t *const p1,
				     const nb_point2D_t *const p2);
static double get_min_dd_in_container(const nb_container_t *const cnt,
				      const nb_point2D_t *const p1,
				      const nb_point2D_t *const p2);
static double get_min_dd_in_bins(const nb_bins2D_t *const bins,
				 const nb_point2D_t *const p1,
				 const nb_point2D_t *const p2);
static bool all_inside_circle(const nb_container_t *const cnt,
			      const double center[2],
			      double radius);

void cunit_nb_geometric_bot_bins2D(void)
{
	CU_pSuite suite = CU_add_suite("nb/geometric_bot/knn/bins2D.c",
				       suite_init, suite_clean);
	CU_add_test(suite, "create()", test_create);
	CU_add_test(suite, "destroy()", test_destroy);
	CU_add_test(suite, "clear()", test_clear);
	CU_add_test(suite, "insert()", test_insert);
	CU_add_test(suite, "delete()", test_delete);
	CU_add_test(suite, "delete_first()", test_delete_first);
	CU_add_test(suite, "get_knn()", test_get_knn);
	CU_add_test(suite, "set_filter()", test_set_filter);
	CU_add_test(suite, "set_filter_data()", test_set_filter_data);
	CU_add_test(suite, "get_candidate_points_to_min_delaunay()",
		    test_get_candidate_points_to_min_delaunay);
	CU_add_test(suite, "get_points_inside_circle()",
		    test_get_points_inside_circle);
	CU_add_test(suite, "are_points_inside_circle()",
		    test_are_points_inside_circle);
	CU_add_test(suite, "get_N_bins()", test_get_N_bins);
	CU_add_test(suite, "get_min_points_x_bin()",
		    test_get_min_points_x_bin);
	CU_add_test(suite, "get_length()", test_get_length);
	CU_add_test(suite, "is_empty()", test_is_empty);
	CU_add_test(suite, "is_not_empty()", test_is_not_empty);
	CU_add_test(suite, "get_size_of_bins()", test_get_size_of_bins);
}

static int suite_init(void)
{
	return 0;
}

static int suite_clean(void)
{
	return 0;
}

static void test_create(void)
{
	nb_bins2D_t *bins = nb_bins2D_create(1.0);
	bool is_ok = (NULL != bins);
	if (is_ok) {
		bool size_is_ok = 
			(fabs(1.0 - nb_bins2D_get_size_of_bins(bins)) < 1e-9);
		bool length_is_ok = (0 == nb_bins2D_get_length(bins));
		is_ok = is_ok && size_is_ok && length_is_ok;
		nb_bins2D_destroy(bins);
	}
	CU_ASSERT(is_ok);
}

static void test_destroy(void)
{
	nb_bins2D_t *bins = nb_bins2D_create(1.0);
	nb_bins2D_destroy(bins);
	CU_ASSERT(true);
}

static void test_clear(void)
{
	int N = 100;
	nb_bins2D_t *bins = get_bins(N);
	nb_bins2D_clear(bins);
	bool is_ok = (0 == nb_bins2D_get_length(bins));
	nb_bins2D_destroy(bins);
	CU_ASSERT(is_ok);
}

static void test_insert(void)
{
	int N = 100;
	nb_bins2D_t *bins = get_bins(N);
	bool is_ok = (N == nb_bins2D_get_length(bins));
	nb_bins2D_destroy(bins);
	CU_ASSERT(is_ok);
}

static void test_delete(void)
{
	int N = 100;
	nb_bins2D_t *bins = get_bins(N);
	nb_point2D_t point;
	point.x[0] = 0.0;
	point.x[1] = 0.0;
	nb_point2D_t *deleted_point = nb_bins2D_delete(bins, &point);
	bool is_ok = false;
	if (NULL != deleted_point) {
		is_ok = (fabs(deleted_point->x[0] - point.x[0]) < 1e-9) &&
			(fabs(deleted_point->x[1] - point.x[1]) < 1e-9);
		nb_point2D_destroy(deleted_point);
	}
	nb_bins2D_destroy(bins);
	CU_ASSERT(is_ok);
}

static void test_delete_first(void)
{
	int N = 100;
	nb_bins2D_t *bins = get_bins(N);
	nb_point2D_t *deleted_point = nb_bins2D_delete_first(bins);
	bool is_ok = false;
	if (NULL != deleted_point) {
		nb_point2D_destroy(deleted_point);
		is_ok = true;
	}
	nb_bins2D_destroy(bins);
	CU_ASSERT(is_ok);
}

static void test_get_knn(void)
{
	int N = 100;
	nb_bins2D_t *bins = get_bins(N);
	nb_point2D_t p;
	p.x[0] = 0.5;
	p.x[1] = 0.5;
	nb_point2D_t* knn[10];
	double knn_dist[10];
	uint32_t k = nb_bins2D_get_knn(bins, &p, 10, knn, knn_dist);
	bool is_ok = (10 == k) && distances_are_ok(10, knn_dist);
	nb_bins2D_destroy(bins);
	CU_ASSERT(is_ok);
}

static void test_set_filter(void)
{
	int N = 100;
	nb_bins2D_t *bins = get_bins(N);
	nb_bins2D_set_filter(bins, filter);
	nb_point2D_t p;
	p.x[0] = 0.5;
	p.x[1] = 0.5;
	nb_point2D_t* knn[10];
	double knn_dist[10];
	uint32_t k = nb_bins2D_get_knn(bins, &p, 10, knn, knn_dist);
	bool is_ok = (10 == k) && distances_are_ok(10, knn_dist) &&
		all_pass_filter(10, &p, knn);
	nb_bins2D_destroy(bins);
	CU_ASSERT(is_ok);
}

static void test_set_filter_data(void)
{
	int N = 100;
	nb_bins2D_t *bins = get_bins(N);
	nb_bins2D_set_filter(bins, filter_with_data);
	int fdata = 5;
	nb_bins2D_set_filter_data(bins, &fdata);
	nb_point2D_t p;
	p.x[0] = 0.5;
	p.x[1] = 0.5;
	nb_point2D_t* knn[10];
	double knn_dist[10];
	uint32_t k = nb_bins2D_get_knn(bins, &p, 10, knn, knn_dist);
	bool is_ok = (10 == k) && distances_are_ok(10, knn_dist) &&
		all_pass_filter(10, &p, knn);
	nb_bins2D_destroy(bins);
	CU_ASSERT(is_ok);
}

static void test_get_candidate_points_to_min_delaunay(void)
{
	int N = 100;
	nb_bins2D_t *bins = get_bins(N);
	nb_point2D_t p1, p2;
	p1.x[0] = -1;
	p1.x[1] = 1;
	p2.x[0] = 1;
	p2.x[1] = -1;

	nb_container_t *cnt = nb_allocate_on_stack(
		nb_container_get_memsize(NB_QUEUE)
	);
	nb_container_init(cnt, NB_QUEUE);
	nb_bins2D_vertices_t vtx_queue;
	vtx_queue.vtx = cnt;
	vtx_queue.add =
		(void (*)(void*, const void *const)) nb_container_insert;
	vtx_queue.is_empty =
		(bool (*)(const void* const)) nb_container_is_empty;

	nb_bins2D_get_candidate_points_to_min_delaunay(bins, &p1, &p2,
						       &vtx_queue);
	bool is_ok = nb_container_is_not_empty(cnt);
	if (is_ok)
		is_ok = is_ok && dont_include_edge_points(cnt,
							  &p1, &p2);
	if (is_ok)
		is_ok = is_ok && all_in_half_space(cnt, &p1, &p2);
	if (is_ok)
		is_ok = is_ok && contains_min_dd(bins, cnt, &p1, &p2);
	nb_container_finish(cnt);
	nb_bins2D_destroy(bins);
	CU_ASSERT(is_ok);
}

static void test_get_points_inside_circle(void)
{
	int N = 100;
	nb_bins2D_t *bins = get_bins(N);
	double center[2] = {0, 0};
	double radius = 10.0;
	nb_container_t *cnt = nb_allocate_on_stack(nb_container_get_memsize(NB_QUEUE));
	nb_container_init(cnt, NB_QUEUE);

	nb_bins2D_get_points_inside_circle(bins, center, radius, cnt);

	bool is_ok = nb_container_is_not_empty(cnt) &&
		all_inside_circle(cnt, center, radius);
	nb_container_finish(cnt);
	nb_bins2D_destroy(bins);
	CU_ASSERT(is_ok);
}

static void test_are_points_inside_circle(void)
{
	int N = 100;
	nb_bins2D_t *bins = get_bins(N);
	double center[2] = {0, 0};
	double radius = 10.0;
	bool is_ok =
		nb_bins2D_are_points_inside_circle(bins, center,
						    radius);
	nb_bins2D_destroy(bins);
	CU_ASSERT(is_ok);
}

static void test_get_N_bins(void)
{
	int N = 100;
	nb_bins2D_t *bins = get_bins(N);
	uint32_t N_bins = nb_bins2D_get_N_bins(bins);
	bool is_ok = (N_bins > 0 && N_bins < N);
	nb_bins2D_destroy(bins);
	CU_ASSERT(is_ok);
}

static void test_get_min_points_x_bin(void)
{
	int N = 100;
	nb_bins2D_t *bins = get_bins(N);
	bool is_ok = (nb_bins2D_get_min_points_x_bin(bins) == 1);
	nb_bins2D_destroy(bins);
	CU_ASSERT(is_ok);
}

static void test_get_length(void)
{
	int N = 100;
	nb_bins2D_t *bins = get_bins(N);
	bool is_ok = (nb_bins2D_get_length(bins) == N);
	nb_bins2D_destroy(bins);
	CU_ASSERT(is_ok);
}

static void test_is_empty(void)
{
	int N = 100;
	nb_bins2D_t *bins = get_bins(N);
	bool is_ok = !nb_bins2D_is_empty(bins);
	nb_bins2D_destroy(bins);
	CU_ASSERT(is_ok);
}

static void test_is_not_empty(void)
{
	int N = 100;
	nb_bins2D_t *bins = get_bins(N);
	bool is_ok = nb_bins2D_is_not_empty(bins);
	nb_bins2D_destroy(bins);
	CU_ASSERT(is_ok);
}

static void test_get_size_of_bins(void)
{
	int N = 100;
	nb_bins2D_t *bins = get_bins(N);
	bool is_ok = (fabs(nb_bins2D_get_size_of_bins(bins) - 1.0) < 1e-9);
	nb_bins2D_destroy(bins);
	CU_ASSERT(is_ok);
}

static nb_bins2D_t* get_bins(int N)
{
	nb_bins2D_t *bins = nb_bins2D_create(1.0);
	for (int i = 0; i < N; i++) {
		nb_point2D_t *point = nb_point2D_create();
		double xsign = (0 == (i/2 + 1) % 2)?(1.0):(-1.0);
		double ysign = (1 == (i/2) % 2)?(1.0):(-1.0);
		point->x[0] = i * 0.5 * xsign;
		point->x[1] = i * 0.5 * ysign;
		nb_bins2D_insert(bins, point);
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

static inline bool filter(const nb_point2D_t *const p_ref,
			  const nb_point2D_t *const p,
			  const void *const data)
{
	return (nb_utils2D_get_dist(p_ref->x, p->x) > 5.0);
}

static inline bool all_pass_filter(int N, nb_point2D_t *p_ref,
				   nb_point2D_t *knn[])
{
	bool is_ok = true;
	for (int i = 0; i < N && is_ok; i++)
		is_ok = is_ok && filter(p_ref, knn[i], NULL);
	return is_ok;
}

static inline bool filter_with_data(const nb_point2D_t *const p_ref,
				    const nb_point2D_t *const p,
				    const void *const data)
{
	return (nb_utils2D_get_dist(p_ref->x, p->x) > *((int*)data));
}

static bool all_in_half_space(const nb_container_t *const cnt,
			      const nb_point2D_t *const p1,
			      const nb_point2D_t *const p2)
{
	bool is_ok = true;
	nb_iterator_t *iter = nb_iterator_create();
	nb_iterator_set_container(iter, cnt);
	while (nb_iterator_has_more(iter)) {
		const nb_point2D_t *p3 = nb_iterator_get_next(iter);
		if (nb_utils2D_orient(p1->x, p2->x, p3->x) < 0.0) {
			is_ok = false;
			break;
		}
	}
	nb_iterator_destroy(iter);
	return is_ok;
}

static bool dont_include_edge_points(const nb_container_t *const cnt,
				     const nb_point2D_t *const p1,
				     const nb_point2D_t *const p2)
{
	bool is_ok = true;
	nb_iterator_t *iter = nb_iterator_create();
	nb_iterator_set_container(iter, cnt);
	while (nb_iterator_has_more(iter)) {
		const nb_point2D_t *p3 = nb_iterator_get_next(iter);
		if (p1 == p3 || p2 == p3) {
			is_ok = false;
			break;
		}
	}
	nb_iterator_destroy(iter);
	return is_ok;
}

static bool contains_min_dd(const nb_bins2D_t *const bins,
			    const nb_container_t *const cnt,
			    const nb_point2D_t *const p1,
			    const nb_point2D_t *const p2)
{
        double min_dd = get_min_dd_in_container(cnt, p1, p2);
	double min_dd_in_bins = get_min_dd_in_bins(bins, p1, p2);
	return min_dd <= min_dd_in_bins;
}

static double get_min_dd_in_container(const nb_container_t *const cnt,
				      const nb_point2D_t *const p1,
				      const nb_point2D_t *const p2)
{
	double min_dd = 1e30;
	nb_iterator_t *iter = nb_iterator_create();
	nb_iterator_set_container(iter, cnt);
	while (nb_iterator_has_more(iter)) {
		const nb_point2D_t *p3 = nb_iterator_get_next(iter);
		double min = nb_utils2D_get_delaunay_dist(p1->x, p2->x,
							   p3->x);
		if (min < min_dd)
			min_dd = min;
	}
	nb_iterator_destroy(iter);
	return min_dd;
}

static double get_min_dd_in_bins(const nb_bins2D_t *const bins,
				 const nb_point2D_t *const p1,
				 const nb_point2D_t *const p2)
{
	double min_dd = 1e30;
	nb_bins2D_iter_t *iter = nb_bins2D_iter_create();
	nb_bins2D_iter_set_bins(iter, bins);
	while (nb_bins2D_iter_has_more(iter)) {
		const nb_point2D_t *p3 = nb_bins2D_iter_get_next(iter);
		if (p1 != p3 && p2 != p3) {
			double orient =
				nb_utils2D_orient(p1->x, p2->x, p3->x);
			if (orient > 0.0) {
				double min = 
					nb_utils2D_get_delaunay_dist(p1->x,
								      p2->x,
								      p3->x);
				if (min < min_dd)
					min_dd = min;
			}
		}
	}
	nb_bins2D_iter_destroy(iter);
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
		const nb_point2D_t *p = nb_iterator_get_next(iter);
		double dist2 = nb_utils2D_get_dist2(center, p->x);
		if (POW2(radius) - dist2 < 0) {
			all_inside = false;
			break;
		}		
	}
	nb_iterator_destroy(iter);
	return all_inside;
}
