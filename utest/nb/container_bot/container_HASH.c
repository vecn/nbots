#define CONTAINER_ID NB_HASH
#define N_ITEMS 200
#include "container_ALL.h"

static int suite_init(void);
static int suite_clean(void);

static void test_realloc_if_max_load(void);
static void test_get_size(void);
static void test_get_N_collisions(void);
static void test_get_collisions(void);
static void test_exist_if_not_with_collision(void);

static uint32_t small_hash_key(const void *const val_ptr);
static bool check_collisions(nb_container_t *collisions);

void cunit_nb_container_bot_HASH(void)
{
	CU_pSuite suite = CU_add_suite("nb/container_bot/container_HASH.c",
				       suite_init, suite_clean);
	container_add_tests(suite);
	CU_add_test(suite, "realloc() if max load",
		    test_realloc_if_max_load);
	CU_add_test(suite, "do('get_size')", test_get_size);
	CU_add_test(suite, "do('get_N_collisions')",
		    test_get_N_collisions);
	CU_add_test(suite, "do('get_collisions')",
		    test_get_collisions);
	CU_add_test(suite, "exist() if not exist and the key collides",
		    test_exist_if_not_with_collision);
}

static int suite_init(void)
{
	return 0;
}

static int suite_clean(void)
{
	return 0;
}

static void test_realloc_if_max_load(void)
{
	int N = 3000;
	nb_container_t *cnt = get_container(N);
	uint32_t length = nb_container_get_length(cnt);
	nb_container_destroy(cnt);
	CU_ASSERT(N == length);  
}

static void test_get_size(void)
{
	nb_container_t *cnt = nb_container_create(NB_HASH);
	int8_t status;
	uint32_t *size_ptr = nb_container_do(cnt, "get_size", NULL, &status);
	uint32_t size = 0;
	if (0 == status) {
		size = *size_ptr;
		free(size_ptr);
	}
	nb_container_destroy(cnt);
	CU_ASSERT(0 < size);
}

static void test_get_N_collisions(void)
{
	nb_container_t *cnt = nb_container_create(NB_HASH);
	nb_container_set_key_generator(cnt, small_hash_key);
	nb_container_set_destroyer(cnt, free);
	for (int i = 0; i < 100; i++) {
		int32_t *val = malloc(sizeof(*val));
		*val = i;
		nb_container_insert(cnt, val);
	}
	int8_t status;
	uint32_t *N_collisions_ptr =
		nb_container_do(cnt, "get_N_collisions", NULL, &status);
	uint32_t N_collisions = 0;
	if (0 == status) {
		N_collisions = *N_collisions_ptr;
		free(N_collisions_ptr);
	}
	CU_ASSERT(3 == N_collisions);
}

static void test_get_collisions(void)
{
	nb_container_t *cnt = nb_container_create(NB_HASH);
	nb_container_set_key_generator(cnt, small_hash_key);
	nb_container_set_destroyer(cnt, free);
	for (int i = 0; i < 100; i++) {
		int32_t *val = malloc(sizeof(*val));
		*val = i;
		nb_container_insert(cnt, val);
	}
	int8_t status;
	nb_container_t *collisions =
		nb_container_do(cnt, "get_collisions", NULL, &status);
	bool is_ok = false;
	if (0 == status) {
		is_ok = check_collisions(collisions);
		nb_container_destroy(collisions);
	}
	CU_ASSERT(is_ok);
}


static void test_exist_if_not_with_collision(void)
{
	int N = N_ITEMS;
	nb_container_t *cnt = get_container(N);
	nb_container_set_comparer(cnt, compare);
	int8_t status;
	uint32_t *size = nb_container_do(cnt, "get_size", NULL, &status);
	int32_t to_find = *size + N_ITEMS; /* Not exist and key collides */
	free(size);
	int32_t *val = nb_container_exist(cnt, &to_find);
	bool is_ok = (NULL == val);
	nb_container_destroy(cnt);
	CU_ASSERT(is_ok);
}

static uint32_t small_hash_key(const void *const val_ptr)
{
	const int32_t *const val = val_ptr;
	return *val % 3;
}

static bool check_collisions(nb_container_t *collisions)
{
	bool is_ok = true;
	while (nb_container_is_not_empty(collisions)) {
		nb_container_t *items = 
			nb_container_delete_first(collisions);
		
		int32_t *val = nb_container_delete_first(items);
		uint32_t key = small_hash_key(val);
		while (nb_container_is_not_empty(items)) {
			val = nb_container_delete_first(items);
			if (small_hash_key(val) != key) {
				is_ok = false;
				break;
			}
		}
		nb_container_destroy(items);
	}
	return is_ok;
}
