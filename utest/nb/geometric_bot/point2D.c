#include <stdbool.h>
#include <stdlib.h>

#include <CUnit/Basic.h>

#include "nb/geometric_bot/point2D.h"

static int suite_init(void);
static int suite_clean(void);

static void test_create(void);
static void test_destroy(void);
static void test_compare(void);

void cunit_nb_geometric_bot_point2D(void)
{
	CU_pSuite suite = CU_add_suite("nb/geometric_bot/point2D.c",
				       suite_init, suite_clean);
	CU_add_test(suite, "create()", test_create);
	CU_add_test(suite, "destroy()", test_destroy);
	CU_add_test(suite, "compare()", test_compare);
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
	nb_point2D_t *p = nb_point2D_create();
	bool is_ok = (0 == p->x[0]) && (0 == p->x[1]) &&
		(NULL == p->attr);
	nb_point2D_destroy(p);
	CU_ASSERT(is_ok);
}

static void test_destroy(void)
{
	nb_point2D_t *p = nb_point2D_create();
	nb_point2D_destroy(p);
	CU_ASSERT(true);
}

static void test_compare(void)
{
	nb_point2D_t p1;
	nb_point2D_t p2;
	p1.x[0] = 1.0 + 1e-16;
	p1.x[1] = 1.0 + 1e-16;
	p2.x[0] = 1.0 - 1e-16;
	p2.x[1] = 1.0 - 1e-16;
	CU_ASSERT(0 == nb_point2D_compare(&p1, &p2));
}
