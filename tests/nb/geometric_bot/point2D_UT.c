#include <stdbool.h>
#include <stdlib.h>

#include "nb/geometric_bot/point2D.h"

#include "test_library.h"
#include "test_add.h"

static bool check_create(void);
static bool check_destroy(void);
static bool check_are_equal(void);

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
	vcn_test_add(tests_ptr, check_are_equal,
		     "Check are_equal()");
}

static bool check_create(void)
{
	vcn_point2D_t *p = vcn_point2D_create();
	bool is_ok = (0 == p->x[0]) && (0 == p->x[1]) &&
		(NULL == p->attr);
	vcn_point2D_destroy(p);
	return is_ok;
}

static bool check_destroy(void)
{
	vcn_point2D_t *p = vcn_point2D_create();
	vcn_point2D_destroy(p);
	return true;
}

static bool check_are_equal(void)
{
	vcn_point2D_t p1;
	vcn_point2D_t p2;
	p1.x[0] = 1.0 + 1e-16;
	p1.x[1] = 1.0 + 1e-16;
	p2.x[0] = 1.0 - 1e-16;
	p2.x[1] = 1.0 - 1e-16;
	return vcn_point2D_are_equal(&p1, &p2);
}
