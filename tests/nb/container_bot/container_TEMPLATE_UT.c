/*
 * This file should not be compiled alone.
 * It must be included after defining the macros
 * "CONTAINER_ID" and "N_ITEMS".
 *
 * Example:
 *  #define CONTAINER_ID NB_CONTAINER_QUEUE
 *  #define N_ITEMS 1000
 *  #include "container_TEMPLATE_UT.c"
 */
#include <stdbool.h>
#include <stdlib.h>
#include <stdint.h>

#include "nb/container_bot/container.h"
#include "nb/statistics_bot.h"

#include "test_library.h"
#include "test_add.h"

static bool check_create(void);
static bool check_clone(void);
static bool check_merge_with_QUEUE(void);
static bool check_merge_with_STACK(void);
static bool check_merge_with_SORTED(void);
static bool check_merge_with_HEAP(void);
static bool check_merge_with_HASH(void);
static bool check_cast_to_QUEUE(void);
static bool check_cast_to_STACK(void);
static bool check_cast_to_SORTED(void);
static bool check_cast_to_HEAP(void);
static bool check_cast_to_HASH(void);
static bool check_cast_to_array(void);
static bool check_destroy(void);
static bool check_clear(void);
static bool check_copy_to_array(void);
static bool check_set_key_generator(void);
static bool check_set_destroyer(void);
static bool check_set_comparer(void);
static bool check_set_cloner(void);
static bool check_insert(void);
static bool check_insert_repeated_item(void);
static bool check_insert_array(void);
static bool check_get_first(void);
static bool check_delete_first(void);
static bool check_exist(void);
static bool check_exist_if_not(void);
static bool check_delete(void);
static bool check_get_length(void);
static bool check_is_empty(void);
static bool check_is_not_empty(void);
static bool check_get_id(void);
static bool check_invalid_do(void);
static bool check_get_dst(void);

static bool check_merge_with(int id);
static bool check_cast_to(int id);
static vcn_container_t* get_container_by_id(int N, int id);
static bool insert_N_int32(vcn_container_t *cnt, int N);
static vcn_container_t* get_container(int N);
static uint32_t keygen(const void *const val);
static bool are_equal(const void *const a, const void *const b);
static void* clone_int32(const void *const val);
static void** get_N_array(vcn_container_t *cnt, int N);
static bool first_is_ok(int32_t *val);

inline int TEMPLATE_get_driver_id(void)
{
	return NB_DRIVER_UNIT_TEST;
}

void TEMPLATE_load_tests(void *tests_ptr)
{
	vcn_test_add(tests_ptr, check_create,
		     "Check create()");
	vcn_test_add(tests_ptr, check_clone,
		     "Check clone()");
	vcn_test_add(tests_ptr, check_merge_with_QUEUE,
		     "Check merge() with QUEUE");
	vcn_test_add(tests_ptr, check_merge_with_STACK,
		     "Check merge() with STACK");
	vcn_test_add(tests_ptr, check_merge_with_SORTED,
		     "Check merge() with SORTED");
	vcn_test_add(tests_ptr, check_merge_with_HEAP,
		     "Check merge() with HEAP");
	vcn_test_add(tests_ptr, check_merge_with_HASH,
		     "Check merge() with HASH");
	vcn_test_add(tests_ptr, check_cast_to_QUEUE,
		     "Check cast() to QUEUE");
	vcn_test_add(tests_ptr, check_cast_to_STACK,
		     "Check cast() to STACK");
	vcn_test_add(tests_ptr, check_cast_to_SORTED,
		     "Check cast() to SORTED");
	vcn_test_add(tests_ptr, check_cast_to_HEAP,
		     "Check cast() to HEAP");
	vcn_test_add(tests_ptr, check_cast_to_HASH,
		     "Check cast() to HASH");
	vcn_test_add(tests_ptr, check_cast_to_array,
		     "Check cast_to_array()");
	vcn_test_add(tests_ptr, check_destroy,
		     "Check destroy()");
	vcn_test_add(tests_ptr, check_clear,
		     "Check clear()");
	vcn_test_add(tests_ptr, check_copy_to_array,
		     "Check copy_to_array()");
	vcn_test_add(tests_ptr, check_set_key_generator,
		     "Check set_key_generator()");
	vcn_test_add(tests_ptr, check_set_destroyer,
		     "Check set_destroyer()");
	vcn_test_add(tests_ptr, check_set_comparer,
		     "Check set_comparer()");
	vcn_test_add(tests_ptr, check_set_cloner,
		     "Check set_cloner()");
	vcn_test_add(tests_ptr, check_insert,
		     "Check insert()");
	vcn_test_add(tests_ptr, check_insert_repeated_item,
		     "Check insert() repeated items");
	vcn_test_add(tests_ptr, check_insert_array,
		     "Check insert_array()");
	vcn_test_add(tests_ptr, check_get_first,
		     "Check get_first()");
	vcn_test_add(tests_ptr, check_delete_first,
		     "Check delete_first()");
	vcn_test_add(tests_ptr, check_exist,
		     "Check exist()");
	vcn_test_add(tests_ptr, check_exist_if_not,
		     "Check exist() if not");
	vcn_test_add(tests_ptr, check_delete,
		     "Check delete()");
	vcn_test_add(tests_ptr, check_get_length,
		     "Check get_length()");
	vcn_test_add(tests_ptr, check_is_empty,
		     "Check is_empty()");
	vcn_test_add(tests_ptr, check_is_not_empty,
		     "Check is_not_empty()");
	vcn_test_add(tests_ptr, check_get_id,
		     "Check get_id()");
	vcn_test_add(tests_ptr, check_invalid_do,
		     "Check do('invalid function')");
	vcn_test_add(tests_ptr, check_get_dst,
		     "Check get_dst()");
}

static bool check_create(void)
{
	vcn_container_t *cnt = vcn_container_create(CONTAINER_ID);
	bool is_ok = (NULL != cnt);
	if (is_ok)
	  vcn_container_destroy(cnt);
	return is_ok;
}

static bool check_clone(void)
{
	int N = N_ITEMS;
	vcn_container_t *cnt = get_container(N);
	vcn_container_set_cloner(cnt, clone_int32);
	vcn_container_t *cloned = vcn_container_clone(cnt);
	uint32_t length = vcn_container_get_length(cloned);
	vcn_container_destroy(cnt);
	bool is_ok = (vcn_container_get_id(cloned) == CONTAINER_ID);
	vcn_container_destroy(cloned);
	return (N == length) && is_ok; 
}

static inline bool check_merge_with_QUEUE(void)
{
	return check_merge_with(NB_CONTAINER_QUEUE);
}

static inline bool check_merge_with_STACK(void)
{
	return check_merge_with(NB_CONTAINER_STACK);
}

static inline bool check_merge_with_SORTED(void)
{
	return check_merge_with(NB_CONTAINER_SORTED);
}

static inline bool check_merge_with_HEAP(void)
{
	return check_merge_with(NB_CONTAINER_HEAP);
}

static inline bool check_merge_with_HASH(void)
{
	return check_merge_with(NB_CONTAINER_HASH);
}

static inline bool check_cast_to_QUEUE(void)
{
	return check_cast_to(NB_CONTAINER_QUEUE);
}

static inline bool check_cast_to_STACK(void)
{
	return check_cast_to(NB_CONTAINER_STACK);
}

static inline bool check_cast_to_SORTED(void)
{
	return check_cast_to(NB_CONTAINER_SORTED);
}

static inline bool check_cast_to_HEAP(void)
{
	return check_cast_to(NB_CONTAINER_HEAP);
}

static inline bool check_cast_to_HASH(void)
{
	return check_cast_to(NB_CONTAINER_HASH);
}

static bool check_cast_to_array(void)
{
	int N = N_ITEMS;
	vcn_container_t *cnt = get_container(N);
	void **array = vcn_container_cast_to_array(cnt);
	for (int i = 0; i < N; i++)
		free(array[i]);
	free(array);
	return true;
}

static bool check_destroy(void)
{
	int N = N_ITEMS;
	vcn_container_t *cnt = get_container(N);
	vcn_container_destroy(cnt);
	return true;
}

static bool check_clear(void)
{
	int N = N_ITEMS;
	vcn_container_t *cnt = get_container(N);
	vcn_container_clear(cnt);
	bool is_empty = vcn_container_is_empty(cnt);
	vcn_container_destroy(cnt);
	return is_empty;
}

static bool check_copy_to_array(void)
{
	int N = N_ITEMS;
	vcn_container_t *cnt = get_container(N);
	void **array = malloc(N * sizeof(*array));
	vcn_container_copy_to_array(cnt, array);
	uint32_t length = vcn_container_get_length(cnt);
	vcn_container_destroy(cnt);
	free(array);
	return (N == length);
}

static bool check_set_key_generator(void)
{
	int N = N_ITEMS;
	vcn_container_t *cnt = get_container(N);
	vcn_container_set_key_generator(cnt, keygen);
	uint32_t length = vcn_container_get_length(cnt);
	vcn_container_destroy(cnt);
	return (N == length);
}

static bool check_set_destroyer(void)
{
	int N = N_ITEMS;
	vcn_container_t *cnt = get_container(N);
	vcn_container_destroy(cnt);
	return true;
}

static bool check_set_comparer(void)
{
	int N = N_ITEMS;
	vcn_container_t *cnt = get_container(N);
	vcn_container_set_comparer(cnt, are_equal);
	int32_t to_find = N_ITEMS;
	int32_t *val = vcn_container_exist(cnt, &to_find);
	bool exist = false;
	if (NULL != val)
		exist = (to_find == *val);
	vcn_container_destroy(cnt);
	return exist;	
}

static bool check_set_cloner(void)
{
	int N = N_ITEMS;
	vcn_container_t *cnt = get_container(N);
	uint32_t length = vcn_container_get_length(cnt);
	vcn_container_set_cloner(cnt, clone_int32);

	vcn_container_t *cloned = vcn_container_clone(cnt);
	vcn_container_destroy(cnt);
	uint32_t cloned_length = vcn_container_get_length(cloned);
	vcn_container_destroy(cloned);

	return (cloned_length == length);
}

static bool check_insert(void)
{
	int N = N_ITEMS;
	vcn_container_t *cnt = get_container(N);
	uint32_t length = vcn_container_get_length(cnt);
	vcn_container_destroy(cnt);
	return (N == length);
}

static bool check_insert_repeated_item(void)
{
	int N = N_ITEMS;
	vcn_container_t *cnt = get_container(N);
	bool success = insert_N_int32(cnt, N);/* Again */
	uint32_t length = vcn_container_get_length(cnt);
	vcn_container_destroy(cnt);
	return (2 * N == length) && success;
}

static bool check_insert_array(void)
{
	int N = N_ITEMS;
	vcn_container_t *cnt = vcn_container_create(CONTAINER_ID);
	void **array = get_N_array(cnt, N);
	vcn_container_insert_array(cnt, N, array);
	uint32_t length = vcn_container_get_length(cnt);
	vcn_container_destroy(cnt);
	free(array);
	return (N == length);
}

static bool check_get_first(void)
{
	int N = N_ITEMS;
	vcn_container_t *cnt = get_container(N);
	int32_t *val = vcn_container_get_first(cnt);
	bool is_ok = first_is_ok(val);
	vcn_container_destroy(cnt);
	return is_ok;
}

static bool check_delete_first(void)
{
	int N = N_ITEMS;
	vcn_container_t *cnt = get_container(N);
	int32_t *val = vcn_container_delete_first(cnt);
	bool is_ok = first_is_ok(val);
	free(val);
	vcn_container_destroy(cnt);
	return is_ok;
}

static bool check_exist(void)
{
	int N = N_ITEMS;
	vcn_container_t *cnt = get_container(N);
	vcn_container_set_comparer(cnt, are_equal);
	int32_t to_find = N_ITEMS;
	int32_t *val = vcn_container_exist(cnt, &to_find);
	bool is_ok = false;
	if (NULL != val)
		is_ok = (to_find == *val);
	vcn_container_destroy(cnt);
	return is_ok;
}

static bool check_exist_if_not(void)
{
	int N = N_ITEMS;
	vcn_container_t *cnt = get_container(N);
	vcn_container_set_comparer(cnt, are_equal);
	int32_t to_find = 2 * N_ITEMS; /* Not exist */
	int32_t *val = vcn_container_exist(cnt, &to_find);
	bool is_ok = (NULL == val);
	vcn_container_destroy(cnt);
	return is_ok;
}

static bool check_delete(void)
{
	int N = N_ITEMS;
	vcn_container_t *cnt = get_container(N);
	vcn_container_set_comparer(cnt, are_equal);
	int32_t to_delete = N_ITEMS;
	int32_t *val = vcn_container_delete(cnt, &to_delete);
	bool is_ok = false;
	if (NULL != val)
		is_ok = (to_delete == *val);
	vcn_container_destroy(cnt);
	return is_ok;
}

static bool check_get_length(void)
{
	int N = N_ITEMS;
	vcn_container_t *cnt = get_container(N);
	uint32_t length = vcn_container_get_length(cnt);
	vcn_container_clear(cnt);
	uint32_t length_zero = vcn_container_get_length(cnt);
	vcn_container_destroy(cnt);
	return (N == length) && (0 == length_zero);
}

static bool check_is_empty(void)
{
	vcn_container_t *cnt = vcn_container_create(CONTAINER_ID);
	bool is_empty = vcn_container_is_empty(cnt);
	vcn_container_destroy(cnt);
	return is_empty;
}

static bool check_is_not_empty(void)
{
	int N = 1;
	vcn_container_t *cnt = get_container(N);
	bool is_not_empty = vcn_container_is_not_empty(cnt);
	vcn_container_destroy(cnt);
	return is_not_empty;
}

static bool check_get_id(void)
{
	vcn_container_t *cnt = vcn_container_create(CONTAINER_ID);
	bool is_ok = (vcn_container_get_id(cnt) == CONTAINER_ID);
	vcn_container_destroy(cnt);
	return is_ok;
}

static bool check_invalid_do(void)
{
	vcn_container_t *cnt = vcn_container_create(CONTAINER_ID);
	int8_t status;
	vcn_container_do(cnt, "invalid", NULL, &status);
	bool is_ok = (0 != status);
	vcn_container_destroy(cnt);
	return is_ok;
}

static bool check_get_dst(void)
{
	vcn_container_t *cnt = vcn_container_create(CONTAINER_ID);
	bool is_not_NULL = (NULL != vcn_container_get_dst(cnt));
	vcn_container_destroy(cnt);
	return is_not_NULL;
}

static bool check_merge_with(int id)
{
	int N = N_ITEMS;
	vcn_container_t *cnt1 = get_container(N);
	vcn_container_t *cnt2 = get_container_by_id(N, id);
	vcn_container_merge(cnt1, cnt2);
	uint32_t length = vcn_container_get_length(cnt1);
	bool is_ok = (vcn_container_get_id(cnt1) == CONTAINER_ID);
	vcn_container_set_destroyer(cnt1, free);
	vcn_container_destroy(cnt1);
	vcn_container_destroy(cnt2);
	return (2 * N == length) && is_ok;
}


static bool check_cast_to(int id)
{
	int N = N_ITEMS;
	vcn_container_t *cnt = get_container(N);
	vcn_container_cast(cnt, id);
	bool is_ok = (vcn_container_get_id(cnt) == id);
	is_ok = is_ok && (N == vcn_container_get_length(cnt));
	vcn_container_destroy(cnt);
	return is_ok;
}

static vcn_container_t* get_container_by_id(int N, int id)
{
	vcn_container_t *cnt = vcn_container_create(id);
      	vcn_container_set_key_generator(cnt, keygen);
      	vcn_container_set_destroyer(cnt, free);
	insert_N_int32(cnt, N);
	return cnt;
}

static bool insert_N_int32(vcn_container_t *cnt, int N)
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
		insertion_ok = insertion_ok && vcn_container_insert(cnt, val);
	}
	return insertion_ok;
}

static inline vcn_container_t* get_container(int N)
{
	return get_container_by_id(N, CONTAINER_ID);
}

static inline uint32_t keygen(const void *const val)
{
	return (uint32_t) *((int32_t*)val);
}

static inline bool are_equal(const void *const a, const void *const b)
{
	int32_t a_val = *((int32_t*)a);
	int32_t b_val = *((int32_t*)b);
	return (a_val == b_val);
}

static inline void* clone_int32(const void *const val)
{
	int32_t *cloned = malloc(sizeof(*cloned));
	*cloned = *((int32_t*)val);
	return cloned;
}

static void** get_N_array(vcn_container_t *cnt, int N)
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
	case NB_CONTAINER_QUEUE:
		is_ok = (1 == *val);
		break;
	case NB_CONTAINER_STACK:
		is_ok = (N_ITEMS == *val);
		break;
	case NB_CONTAINER_SORTED:
		is_ok = (1 == *val);
		break;
	case NB_CONTAINER_HEAP:
		is_ok = (1 == *val);
		break;
	case NB_CONTAINER_HASH:
		is_ok = (NULL != val);
		break;
	default:
		is_ok = false;
	}
	return is_ok;
}