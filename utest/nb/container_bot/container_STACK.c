#define CONTAINER_ID NB_STACK
#define N_ITEMS 200

#include "container_ALL.h"

static int suite_init(void);
static int suite_clean(void);

static void test_insert_first(void);

void cunit_nb_container_bot_STACK(void)
{
	CU_pSuite suite = CU_add_suite("nb/container_bot/container_STACK.c",
				       suite_init, suite_clean);
	container_add_tests(suite);
}


static int suite_init(void)
{
	return 0;
}

static int suite_clean(void)
{
	return 0;
}
