/*
 * This file should not be compiled alone.
 * It must be included after defining the macros
 * "CONTAINER_ID" and "N_ITEMS".
 *
 * Example:
 *  #define CONTAINER_ID NB_QUEUE
 *  #define N_ITEMS 1000
 *  #include "container_all.h"
 */
#include <stdbool.h>
#include <stdlib.h>
#include <stdint.h>

#include <CUnit/Basic.h>

#include "nb/container_bot/container.h"
#include "nb/statistics_bot.h"

static void test_create(void);
static void test_clone(void);
static void test_destroy(void);
static void test_clear(void);
static void test_merge_with_QUEUE(void);
static void test_merge_with_STACK(void);
static void test_merge_with_SORTED(void);
static void test_merge_with_HEAP(void);
static void test_merge_with_HASH(void);
static void test_cast_to_array(void);
static void test_copy_to_array(void);
static void test_set_key_generator(void);
static void test_set_destroyer(void);
static void test_set_comparer(void);
static void test_set_cloner(void);
static void test_insert(void);
static void test_insert_repeated_item(void);
static void test_insert_array(void);
static void test_get_first(void);
static void test_delete_first(void);
static void test_exist(void);
static void test_exist_if_not(void);
static void test_delete(void);
static void test_get_length(void);
static void test_is_empty(void);
static void test_is_not_empty(void);
static void test_get_type(void);
static void test_invalid_do(void);

static void test_merge_with(nb_container_type type);
static nb_container_t* get_container_by_type(int N, nb_container_type type);
static bool insert_N_int32(nb_container_t *cnt, int N);
static nb_container_t* get_container(int N);
static uint32_t keygen(const void *const val);
static int8_t compare(const void *const a, const void *const b);
static void* clone_int32(const void *const val);
static void** get_N_array(nb_container_t *cnt, int N);
static bool first_is_ok(int32_t *val);

static void container_add_tests(CU_pSuite suite)
{
	CU_add_test(suite, "create()", test_create);
	CU_add_test(suite, "clone()", test_clone);
	CU_add_test(suite, "destroy()", test_destroy);
	CU_add_test(suite, "clear()", test_clear);
	CU_add_test(suite, "merge() with QUEUE", test_merge_with_QUEUE);
	CU_add_test(suite, "merge() with STACK", test_merge_with_STACK);
	CU_add_test(suite, "merge() with SORTED", test_merge_with_SORTED);
	CU_add_test(suite, "merge() with HEAP", test_merge_with_HEAP);
	CU_add_test(suite, "merge() with HASH", test_merge_with_HASH);
	CU_add_test(suite, "cast_to_array()", test_cast_to_array);
	CU_add_test(suite, "copy_to_array()", test_copy_to_array);
	CU_add_test(suite, "set_key_generator()", test_set_key_generator);
	CU_add_test(suite, "set_destroyer()", test_set_destroyer);
	CU_add_test(suite, "set_comparer()", test_set_comparer);
	CU_add_test(suite, "set_cloner()", test_set_cloner);
	CU_add_test(suite, "insert()", test_insert);
	CU_add_test(suite, "insert() repeated items",
		    test_insert_repeated_item);
	CU_add_test(suite, "insert_array()", test_insert_array);
	CU_add_test(suite, "get_first()", test_get_first);
	CU_add_test(suite, "delete_first()", test_delete_first);
	CU_add_test(suite, "exist()", test_exist);
	CU_add_test(suite, "exist() if not", test_exist_if_not);
	CU_add_test(suite, "delete()", test_delete);
	CU_add_test(suite, "get_length()", test_get_length);
	CU_add_test(suite, "is_empty()", test_is_empty);
	CU_add_test(suite, "is_not_empty()", test_is_not_empty);
	CU_add_test(suite, "get_type()", test_get_type);
	CU_add_test(suite, "do('invalid function')", test_invalid_do);
}

static void test_create(void)
{
	nb_container_t *cnt = nb_container_create(CONTAINER_ID);
	bool is_ok = (NULL != cnt);
	if (is_ok)
	  nb_container_destroy(cnt);
	CU_ASSERT(is_ok);
}

static void test_clone(void)
{
	int N = N_ITEMS;
	nb_container_t *cnt = get_container(N);
	nb_container_set_cloner(cnt, clone_int32);

	nb_container_t *cloned = nb_container_clone(cnt);
	uint32_t length = nb_container_get_length(cloned);
	nb_container_destroy(cnt);
	bool is_ok = (nb_container_get_type(cloned) == CONTAINER_ID);
	nb_container_destroy(cloned);
	CU_ASSERT(N == length);
	CU_ASSERT(is_ok);
}

static void test_destroy(void)
{
	int N = N_ITEMS;
	nb_container_t *cnt = get_container(N);
	nb_container_destroy(cnt);
	CU_ASSERT(true);
}

static void test_clear(void)
{
	int N = N_ITEMS;
	nb_container_t *cnt = get_container(N);
	nb_container_clear(cnt);
	bool is_empty = nb_container_is_empty(cnt);
	nb_container_destroy(cnt);
	CU_ASSERT(is_empty);
}

static void test_merge_with_QUEUE(void)
{
	test_merge_with(NB_QUEUE);
}

static void test_merge_with_STACK(void)
{
	test_merge_with(NB_STACK);
}

static void test_merge_with_SORTED(void)
{
	test_merge_with(NB_SORTED);
}

static void test_merge_with_HEAP(void)
{
	test_merge_with(NB_HEAP);
}

static void test_merge_with_HASH(void)
{
	test_merge_with(NB_HASH);
}

static void test_cast_to_array(void)
{
	int N = N_ITEMS;
	nb_container_t *cnt = get_container(N);
	void **array = nb_container_cast_to_array(cnt);
	for (int i = 0; i < N; i++)
		free(array[i]);
	free(array);
	CU_ASSERT(true);
}

static void test_copy_to_array(void)
{
	int N = N_ITEMS;
	nb_container_t *cnt = get_container(N);
	void **array = malloc(N * sizeof(*array));
	nb_container_copy_to_array(cnt, array);
	uint32_t length = nb_container_get_length(cnt);
	nb_container_destroy(cnt);
	free(array);
	CU_ASSERT(N == length);
}

static void test_set_key_generator(void)
{
	int N = N_ITEMS;
	nb_container_t *cnt = get_container(N);
	nb_container_set_key_generator(cnt, keygen);
	uint32_t length = nb_container_get_length(cnt);
	nb_container_destroy(cnt);
	CU_ASSERT(N == length);
}

static void test_set_destroyer(void)
{
	int N = N_ITEMS;
	nb_container_t *cnt = get_container(N);
	nb_container_destroy(cnt);
	CU_ASSERT(true);
}

static void test_set_comparer(void)
{
	int N = N_ITEMS;
	nb_container_t *cnt = get_container(N);
	nb_container_set_comparer(cnt, compare);
	int32_t to_find = N_ITEMS;
	int32_t *val = nb_container_exist(cnt, &to_find);
	bool exist = false;
	if (NULL != val)
		exist = (to_find == *val);
	nb_container_destroy(cnt);
	CU_ASSERT(exist);	
}

static void test_set_cloner(void)
{
	int N = N_ITEMS;
	nb_container_t *cnt = get_container(N);
	nb_container_set_cloner(cnt, clone_int32);

	nb_container_t *cloned = nb_container_clone(cnt);
	nb_container_destroy(cnt);
	uint32_t length = nb_container_get_length(cloned);
	nb_container_destroy(cloned);
	CU_ASSERT(N == length);
}

static void test_insert(void)
{
	int N = N_ITEMS;
	nb_container_t *cnt = get_container(N);
	uint32_t length = nb_container_get_length(cnt);
	nb_container_destroy(cnt);
	CU_ASSERT(N == length);
}

static void test_insert_repeated_item(void)
{
	int N = N_ITEMS;
	nb_container_t *cnt = get_container(N);
	bool success = insert_N_int32(cnt, N);/* Again */
	uint32_t length = nb_container_get_length(cnt);
	nb_container_destroy(cnt);
	CU_ASSERT(2 * N == length);
	CU_ASSERT(success);
}

static void test_insert_array(void)
{
	int N = N_ITEMS;
	nb_container_t *cnt = nb_container_create(CONTAINER_ID);
	void **array = get_N_array(cnt, N);
	nb_container_insert_array(cnt, N, array);
	uint32_t length = nb_container_get_length(cnt);
	nb_container_destroy(cnt);
	free(array);
	CU_ASSERT(N == length);
}

static void test_get_first(void)
{
	int N = N_ITEMS;
	nb_container_t *cnt = get_container(N);
	int32_t *val = nb_container_get_first(cnt);
	bool is_ok = first_is_ok(val);
	nb_container_destroy(cnt);
	CU_ASSERT(is_ok);
}

static void test_delete_first(void)
{
	int N = N_ITEMS;
	nb_container_t *cnt = get_container(N);
	int32_t *val = nb_container_delete_first(cnt);
	bool is_ok = first_is_ok(val);
	free(val);
	nb_container_destroy(cnt);
	CU_ASSERT(is_ok);
}

static void test_exist(void)
{
	int N = N_ITEMS;
	nb_container_t *cnt = get_container(N);
	int32_t to_find = N_ITEMS;
	int32_t *val = nb_container_exist(cnt, &to_find);
	bool is_ok = false;
	if (NULL != val)
		is_ok = (to_find == *val);
	nb_container_destroy(cnt);
	CU_ASSERT(is_ok);
}

static void test_exist_if_not(void)
{
	int N = N_ITEMS;
	nb_container_t *cnt = get_container(N);
	int32_t to_find = 2 * N_ITEMS; /* Not exist */
	int32_t *val = nb_container_exist(cnt, &to_find);
	bool is_ok = (NULL == val);
	nb_container_destroy(cnt);
	CU_ASSERT(is_ok);
}

static void test_delete(void)
{
	int N = N_ITEMS;
	nb_container_t *cnt = get_container(N);
	int32_t to_delete = N_ITEMS;
	int32_t *val = nb_container_delete(cnt, &to_delete);
	bool is_ok = false;
	if (NULL != val)
		is_ok = (to_delete == *val);
	nb_container_destroy(cnt);
	CU_ASSERT(is_ok);
}

static void test_get_length(void)
{
	int N = N_ITEMS;
	nb_container_t *cnt = get_container(N);
	uint32_t length = nb_container_get_length(cnt);
	nb_container_clear(cnt);
	uint32_t length_zero = nb_container_get_length(cnt);
	nb_container_destroy(cnt);
	CU_ASSERT(N == length);
	CU_ASSERT(0 == length_zero);
}

static void test_is_empty(void)
{
	nb_container_t *cnt = nb_container_create(CONTAINER_ID);
	bool is_empty = nb_container_is_empty(cnt);
	nb_container_destroy(cnt);
	CU_ASSERT(is_empty);
}

static void test_is_not_empty(void)
{
	int N = 1;
	nb_container_t *cnt = get_container(N);
	bool is_not_empty = nb_container_is_not_empty(cnt);
	nb_container_destroy(cnt);
	CU_ASSERT(is_not_empty);
}

static void test_get_type(void)
{
	nb_container_t *cnt = nb_container_create(CONTAINER_ID);
	bool is_ok = (nb_container_get_type(cnt) == CONTAINER_ID);
	nb_container_destroy(cnt);
	CU_ASSERT(is_ok);
}

static void test_invalid_do(void)
{
	nb_container_t *cnt = nb_container_create(CONTAINER_ID);
	int8_t status;
	nb_container_do(cnt, "invalid", NULL, &status);
	bool is_ok = (0 != status);
	nb_container_destroy(cnt);
	CU_ASSERT(is_ok);
}

static void test_merge_with(nb_container_type type)
{
	int N = N_ITEMS;
	nb_container_t *cnt1 = get_container(N);
	nb_container_t *cnt2 = get_container_by_type(N, type);
	nb_container_merge(cnt1, cnt2);
	uint32_t length = nb_container_get_length(cnt1);
	bool is_ok = (nb_container_get_type(cnt1) == CONTAINER_ID);
	nb_container_set_destroyer(cnt1, free);
	nb_container_destroy(cnt1);
	nb_container_destroy(cnt2);
	CU_ASSERT(2 * N == length);
	CU_ASSERT(is_ok);
}

static nb_container_t* get_container_by_type(int N, nb_container_type type)
{
	nb_container_t *cnt = nb_container_create(type);
      	nb_container_set_key_generator(cnt, keygen);
	nb_container_set_comparer(cnt, compare);
      	nb_container_set_destroyer(cnt, free);

	insert_N_int32(cnt, N);
	return cnt;
}

static bool insert_N_int32(nb_container_t *cnt, int N)
{
	uint32_t seed = 1;
	bool insertion_ok = true;
	for (int i = 0; i < N; i++) {
		int32_t *val = malloc(sizeof(*val));
		if (i == N-1) {
			*val = N_ITEMS;
		} else {
			*val = seed;
			seed = vcn_statistics_lcg(seed);
		}
		insertion_ok = insertion_ok && nb_container_insert(cnt, val);
	}
	return insertion_ok;
}

static inline nb_container_t* get_container(int N)
{
	return get_container_by_type(N, CONTAINER_ID);
}

static inline uint32_t keygen(const void *const val)
{
	return (uint32_t) *((int32_t*)val);
}

static int8_t compare(const void *const a, const void *const b)
{
	int32_t a_val = *((int32_t*)a);
	int32_t b_val = *((int32_t*)b);
	int8_t out;
	if (a_val < b_val)
		out = -1;
	else if (a_val > b_val)
		out = 1;
	else
		out = 0;
	return out;
}

static inline void* clone_int32(const void *const val)
{
	int32_t *cloned = malloc(sizeof(*cloned));
	*cloned = *((int32_t*)val);
	return cloned;
}

static void** get_N_array(nb_container_t *cnt, int N)
{
	void **array = malloc(N * sizeof(*array));
	for (int i = 0; i < N; i++) {
		int32_t *val = malloc(sizeof(*val));
		*val = (N - i);
		array[i] = val;
	}
	return array;
}

static bool first_is_ok(int32_t *val)
{
	bool is_ok;
	switch (CONTAINER_ID) {
	case NB_QUEUE:
		is_ok = (1 == *val);
		break;
	case NB_STACK:
		is_ok = (N_ITEMS == *val);
		break;
	case NB_SORTED:
		is_ok = (1 == *val);
		break;
	case NB_HEAP:
		is_ok = (1 == *val);
		break;
	case NB_HASH:
		is_ok = (NULL != val);
		break;
	default:
		is_ok = false;
	}
	return is_ok;
}
