#include <stdbool.h>
#include <stdlib.h>

#include "nb/geometric_bot/point2D.h"
#include "nb/geometric_bot/knn/bins2D.h"
#include "nb/geometric_bot/knn/bins2D_iterator.h"

#include "test_library.h"
#include "test_add.h"

static bool check_create(void);
static bool check_destroy(void);
static bool check_set_bins(void);
static bool check_has_more(void);
static bool check_get_next(void);
static bool check_restart(void);

static vcn_bins2D_t* get_bins(int N);

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
	vcn_test_add(tests_ptr, check_set_bins,
		     "Check set_bins()");
	vcn_test_add(tests_ptr, check_has_more,
		     "Check has_more()");
	vcn_test_add(tests_ptr, check_get_next,
		     "Check get_next()");
	vcn_test_add(tests_ptr, check_restart,
		     "Check restart()");
}

static bool check_create(void)
{
	vcn_bins2D_iter_t *iter = vcn_bins2D_iter_create();
	bool is_ok = false;
	if (NULL != iter) {
		is_ok = true;
		vcn_bins2D_iter_destroy(iter);
	}
	return is_ok;
}

static bool check_destroy(void)
{
	vcn_bins2D_iter_t *iter = vcn_bins2D_iter_create();
	vcn_bins2D_iter_destroy(iter);
	return true;
}

static bool check_set_bins(void)
{
	int N = 100;
	vcn_bins2D_iter_t *iter = vcn_bins2D_iter_create();
	vcn_bins2D_t *bins = get_bins(N);
	vcn_bins2D_iter_set_bins(iter, bins);
	bool is_ok = vcn_bins2D_iter_has_more(iter);
	vcn_bins2D_iter_destroy(iter);
	vcn_bins2D_destroy(bins);
	return is_ok;
}

static bool check_has_more(void)
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
	return (counter == N);
}

static bool check_get_next(void)
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
	return is_ok;
}

static bool check_restart(void)
{
	int N = 100;
	vcn_bins2D_iter_t *iter = vcn_bins2D_iter_create();
	vcn_bins2D_t *bins = get_bins(N);
	vcn_bins2D_iter_set_bins(iter, bins);
	const vcn_point2D_t *first_point = vcn_bins2D_iter_get_next(iter);
	bool is_ok = (NULL != first_point);
	vcn_bins2D_iter_restart(iter);
	const vcn_point2D_t *first_point_after_restart = vcn_bins2D_iter_get_next(iter);
	is_ok = is_ok && (first_point_after_restart == first_point);
	vcn_bins2D_iter_destroy(iter);
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
