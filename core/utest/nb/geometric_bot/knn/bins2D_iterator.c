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

static vcn_bins2D_t* get_bins(int N);

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
	vcn_bins2D_iter_t *iter = vcn_bins2D_iter_create();
	bool is_ok = false;
	if (NULL != iter) {
		is_ok = true;
		vcn_bins2D_iter_destroy(iter);
	}
	CU_ASSERT(is_ok);
}

static void test_destroy(void)
{
	vcn_bins2D_iter_t *iter = vcn_bins2D_iter_create();
	vcn_bins2D_iter_destroy(iter);
	CU_ASSERT(true);
}

static void test_set_bins(void)
{
	int N = 100;
	vcn_bins2D_iter_t *iter = vcn_bins2D_iter_create();
	vcn_bins2D_t *bins = get_bins(N);
	vcn_bins2D_iter_set_bins(iter, bins);
	bool is_ok = vcn_bins2D_iter_has_more(iter);
	vcn_bins2D_iter_destroy(iter);
	vcn_bins2D_destroy(bins);
	CU_ASSERT(is_ok);
}

static void test_has_more(void)
{
	int N = 100;
	vcn_bins2D_iter_t *iter = vcn_bins2D_iter_create();
	vcn_bins2D_t *bins = get_bins(N);
	vcn_bins2D_iter_set_bins(iter, bins);
	int counter = 0;
	while (vcn_bins2D_iter_has_more(iter)) {
		vcn_bins2D_iter_get_next(iter);
		counter += 1;
	}
	vcn_bins2D_iter_destroy(iter);
	vcn_bins2D_destroy(bins);
	CU_ASSERT(counter == N);
}

static void test_get_next(void)
{
	int N = 100;
	vcn_bins2D_iter_t *iter = vcn_bins2D_iter_create();
	vcn_bins2D_t *bins = get_bins(N);
	vcn_bins2D_iter_set_bins(iter, bins);
	bool is_ok = true;
	while (vcn_bins2D_iter_has_more(iter)) {
		const vcn_point2D_t *point = vcn_bins2D_iter_get_next(iter);
		if (NULL == point) {
			is_ok = false;
			break;
		}
	}
	vcn_bins2D_iter_destroy(iter);
	vcn_bins2D_destroy(bins);
	CU_ASSERT(is_ok);
}

static void test_restart(void)
{
	int N = 100;
	vcn_bins2D_iter_t *iter = vcn_bins2D_iter_create();
	vcn_bins2D_t *bins = get_bins(N);
	vcn_bins2D_iter_set_bins(iter, bins);
	const vcn_point2D_t *first_point =
		vcn_bins2D_iter_get_next(iter);
	bool is_ok = (NULL != first_point);
	vcn_bins2D_iter_restart(iter);
	const vcn_point2D_t *first_point_after_restart =
		vcn_bins2D_iter_get_next(iter);
	is_ok = is_ok && (first_point_after_restart == first_point);
	vcn_bins2D_iter_destroy(iter);
	vcn_bins2D_destroy(bins);
	CU_ASSERT(is_ok);
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
