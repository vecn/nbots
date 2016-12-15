/*
 * This file should not be compiled alone.
 * It must be included after defining the macros
 * "CONTAINER_ID" and "N_ITEMS".
 *
 * Example:
 *  #define CONTAINER_ID NB_QUEUE
 *  #define N_ITEMS 1000
 *  #include "iterator_ALL.c"
 */
#include <stdbool.h>
#include <stdlib.h>
#include <stdint.h>

#include "cunit/Basic.h"

#include "nb/container_bot/container.h"
#include "nb/container_bot/iterator.h"

static void test_create(void);
static void test_clone(void);
static void test_destroy(void);
static void test_set_container(void);
static void test_restart(void);
static void test_get_next(void);
static void test_has_more(void);
static void test_has_more_with_1_item(void);

static nb_container_t* get_container(int N);
static void insert_N_int32(nb_container_t *cnt, int N);
static uint32_t keygen(const void *const val);
static int8_t compare(const void *const a, const void *const b);

static void iterator_add_tests(CU_pSuite suite)
{
	CU_add_test(suite, "create()", test_create);
	CU_add_test(suite, "clone()", test_clone);
	CU_add_test(suite, "destroy()", test_destroy);
	CU_add_test(suite, "set_container()", test_set_container);
	CU_add_test(suite, "restart()", test_restart);
	CU_add_test(suite, "get_next()", test_get_next);
	CU_add_test(suite, "has_more()", test_has_more);
	CU_add_test(suite, "has_more() with one item",
		    test_has_more_with_1_item);
}

static void test_create(void)
{
	nb_iterator_t *iter = nb_iterator_create();
	bool is_ok = (NULL != iter);
	if (is_ok)
		nb_iterator_destroy(iter);
	CU_ASSERT(is_ok);
}

static void test_clone(void)
{
	nb_container_t *cnt = get_container(N_ITEMS);
	nb_iterator_t *iter = nb_iterator_create();
	nb_iterator_set_container(iter, cnt);
	nb_iterator_t *cloned = nb_iterator_clone(iter);
	bool is_ok = nb_iterator_has_more(cloned);
	nb_iterator_destroy(cloned);
	nb_iterator_destroy(iter);
	nb_container_destroy(cnt);
	CU_ASSERT(is_ok);
}

static void test_destroy(void)
{
	nb_iterator_t *iter = nb_iterator_create();
	nb_iterator_destroy(iter);
	CU_ASSERT(true);
}

static void test_set_container(void)
{
	nb_container_t *cnt = get_container(N_ITEMS);
	nb_iterator_t *iter = nb_iterator_create();
	nb_iterator_set_container(iter, cnt);
	bool is_ok = nb_iterator_has_more(iter);
	nb_iterator_destroy(iter);
	nb_container_destroy(cnt);
	CU_ASSERT(is_ok);
}

static void test_restart(void)
{
	nb_container_t *cnt = get_container(N_ITEMS);
	nb_iterator_t *iter = nb_iterator_create();
	nb_iterator_set_container(iter, cnt);
	while (nb_iterator_has_more(iter))
		nb_iterator_get_next(iter);
	nb_iterator_restart(iter);
	bool is_ok = nb_iterator_has_more(iter);
	nb_iterator_destroy(iter);
	nb_container_destroy(cnt);
	CU_ASSERT(is_ok);
}

static void test_get_next(void)
{
	nb_container_t *cnt = get_container(N_ITEMS);
	nb_iterator_t *iter = nb_iterator_create();
	nb_iterator_set_container(iter, cnt);
	bool is_ok = true;
	while (nb_iterator_has_more(iter)) {
		const int32_t *val = nb_iterator_get_next(iter);
		is_ok = is_ok && (NULL != val);
	}
	nb_iterator_destroy(iter);
	nb_container_destroy(cnt);
	CU_ASSERT(is_ok);
}

static void test_has_more(void)
{
	nb_container_t *cnt = get_container(N_ITEMS);
	nb_iterator_t *iter = nb_iterator_create();
	nb_iterator_set_container(iter, cnt);
	int32_t counter = 0;
	while (nb_iterator_has_more(iter)) {
		nb_iterator_get_next(iter);
		counter += 1;
	}
	nb_iterator_destroy(iter);
	nb_container_destroy(cnt);
	CU_ASSERT(N_ITEMS == counter);
}

static void test_has_more_with_1_item(void)
{
	nb_container_t *cnt = get_container(1);
	nb_iterator_t *iter = nb_iterator_create();
	nb_iterator_set_container(iter, cnt);
	int32_t counter = 0;
	while (nb_iterator_has_more(iter)) {
		nb_iterator_get_next(iter);
		counter += 1;
	}
	nb_iterator_destroy(iter);
	nb_container_destroy(cnt);
	CU_ASSERT(1 == counter);
}

static inline nb_container_t* get_container(int N)
{
	nb_container_t *cnt = nb_container_create(CONTAINER_ID);
      	nb_container_set_key_generator(cnt, keygen);
      	nb_container_set_comparer(cnt, compare);
      	nb_container_set_destroyer(cnt, free);
	insert_N_int32(cnt, N);
	return cnt;
}

static void insert_N_int32(nb_container_t *cnt, int N)
{
	for (int32_t i = 0; i < N; i++) {
		int32_t *val = malloc(sizeof(*val));
		*val = i;
		nb_container_insert(cnt, val);
	}
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
	if (a_val - b_val < 0)
		out = -1;
	else if (a_val - b_val > 0)
		out = 1;
	else
		out = 0;
	return out;
}
