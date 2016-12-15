#define CONTAINER_ID NB_QUEUE
#define N_ITEMS 200

#include "container_ALL.h"

static int suite_init(void);
static int suite_clean(void);

static void test_insert_first(void);

void cunit_nb_container_bot_QUEUE(void)
{
	CU_pSuite suite = CU_add_suite("nb/container_bot/container_QUEUE.c",
				       suite_init, suite_clean);
	container_add_tests(suite);
	CU_add_test(suite, "do('insert_first')", test_insert_first);
}

static int suite_init(void)
{
	return 0;
}

static int suite_clean(void)
{
	return 0;
}

static void test_insert_first(void)
{
	int N = N_ITEMS;
	nb_container_t *cnt = nb_container_create(CONTAINER_ID);
	bool is_ok = true;
	for (int i = 0; i < N; i++) {
		int *value = malloc(sizeof(*value));
		*value = i;
		int8_t status;
		nb_container_do(cnt, "insert_first", value, &status);
		is_ok = is_ok && (0 == status);
	}
	int32_t length = nb_container_get_length(cnt);
	int *first_ptr = nb_container_get_first(cnt);
	int first_val = *first_ptr;
	nb_container_set_destroyer(cnt, free);
	nb_container_destroy(cnt);

	CU_ASSERT(is_ok);
	CU_ASSERT(N == length);
	CU_ASSERT((N - 1) == first_val);
}
