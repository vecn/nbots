#include <stdbool.h>
#include <stdlib.h>

#include <CUnit/Basic.h>

#include "nb/geometric_bot/point2D.h"
#include "nb/geometric_bot/knn/bins2D.h"
#include "nb/geometric_bot/knn/bins2D_iterator.h"

static int suite_init(void);
static int suite_clean(void);

static void test_create(void);
static void test_destroy(void);
static void test_set_bins(void);
static void test_has_more(void);
static void test_get_next(void);
static void test_restart(void);

static nb_bins2D_t* get_bins(int N);

void cunit_nb_geometric_bot_bins2D_iterator(void)
{
	CU_pSuite suite = 
	  CU_add_suite("nb/geometric_bot/knn/bins2D_iterator.c",
		       suite_init, suite_clean);
	CU_add_test(suite, "create()", test_create);
	CU_add_test(suite, "destroy()", test_destroy);
	CU_add_test(suite, "set_bins()", test_set_bins);
	CU_add_test(suite, "has_more()", test_has_more);
	CU_add_test(suite, "get_next()", test_get_next);
	CU_add_test(suite, "restart()", test_restart);
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
	nb_bins2D_iter_t *iter = nb_bins2D_iter_create();
	bool is_ok = false;
	if (NULL != iter) {
		is_ok = true;
		nb_bins2D_iter_destroy(iter);
	}
	CU_ASSERT(is_ok);
}

static void test_destroy(void)
{
	nb_bins2D_iter_t *iter = nb_bins2D_iter_create();
	nb_bins2D_iter_destroy(iter);
	CU_ASSERT(true);
}

static void test_set_bins(void)
{
	int N = 100;
	nb_bins2D_iter_t *iter = nb_bins2D_iter_create();
	nb_bins2D_t *bins = get_bins(N);
	nb_bins2D_iter_set_bins(iter, bins);
	bool is_ok = nb_bins2D_iter_has_more(iter);
	nb_bins2D_iter_destroy(iter);
	nb_bins2D_destroy(bins);
	CU_ASSERT(is_ok);
}

static void test_has_more(void)
{
	int N = 100;
	nb_bins2D_iter_t *iter = nb_bins2D_iter_create();
	nb_bins2D_t *bins = get_bins(N);
	nb_bins2D_iter_set_bins(iter, bins);
	int counter = 0;
	while (nb_bins2D_iter_has_more(iter)) {
		nb_bins2D_iter_get_next(iter);
		counter += 1;
	}
	nb_bins2D_iter_destroy(iter);
	nb_bins2D_destroy(bins);
	CU_ASSERT(counter == N);
}

static void test_get_next(void)
{
	int N = 100;
	nb_bins2D_iter_t *iter = nb_bins2D_iter_create();
	nb_bins2D_t *bins = get_bins(N);
	nb_bins2D_iter_set_bins(iter, bins);
	bool is_ok = true;
	while (nb_bins2D_iter_has_more(iter)) {
		const nb_point2D_t *point = nb_bins2D_iter_get_next(iter);
		if (NULL == point) {
			is_ok = false;
			break;
		}
	}
	nb_bins2D_iter_destroy(iter);
	nb_bins2D_destroy(bins);
	CU_ASSERT(is_ok);
}

static void test_restart(void)
{
	int N = 100;
	nb_bins2D_iter_t *iter = nb_bins2D_iter_create();
	nb_bins2D_t *bins = get_bins(N);
	nb_bins2D_iter_set_bins(iter, bins);
	const nb_point2D_t *first_point =
		nb_bins2D_iter_get_next(iter);
	bool is_ok = (NULL != first_point);
	nb_bins2D_iter_restart(iter);
	const nb_point2D_t *first_point_after_restart =
		nb_bins2D_iter_get_next(iter);
	is_ok = is_ok && (first_point_after_restart == first_point);
	nb_bins2D_iter_destroy(iter);
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
