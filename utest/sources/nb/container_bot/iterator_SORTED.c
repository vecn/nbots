#define CONTAINER_ID NB_SORTED
#define N_ITEMS 5000
#include "iterator_ALL.h"

static int suite_init(void);
static int suite_clean(void);

void cunit_nb_container_bot_iterator_SORTED(void)
{
	CU_pSuite suite = CU_add_suite("nb/container_bot/iterator_SORTED.c",
				       suite_init, suite_clean);
	iterator_add_tests(suite);
}

static int suite_init(void)
{
	return 0;
}

static int suite_clean(void)
{
	return 0;
}
