#define CONTAINER_ID NB_SORTED
#define N_ITEMS 200

#include "container_ALL.h"

static int suite_init(void);
static int suite_clean(void);

static void test_delete_last(void);

void cunit_nb_container_bot_SORTED(void)
{
	CU_pSuite suite = CU_add_suite("nb/container_bot/container_SORTED.c",
				       suite_init, suite_clean);
	container_add_tests(suite);
	CU_add_test(suite, "do('delete_last')", test_delete_last);
}

static int suite_init(void)
{
	return 0;
}

static int suite_clean(void)
{
	return 0;
}

static void test_delete_last(void)
{
	int N = N_ITEMS;
	nb_container_t *cnt = nb_container_create(CONTAINER_ID);
	nb_container_set_comparer(cnt, compare);
	nb_container_set_destroyer(cnt, free);
	for (int i = 0; i < N; i++) {
		int *value = malloc(sizeof(*value));
		*value = i;
		nb_container_insert(cnt, value);
	}
	int8_t status;
	int32_t *val = nb_container_do(cnt, "delete_last", NULL, &status);
	bool is_ok = ((N - 1) == *val);
	free(val);
	bool length_ok = ((N - 1) == nb_container_get_length(cnt));
	nb_container_destroy(cnt);
	CU_ASSERT(is_ok);
	CU_ASSERT(length_ok);
}
