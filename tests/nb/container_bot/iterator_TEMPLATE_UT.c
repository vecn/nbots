/*
 * This file should not be compiled alone.
 * It must be included after defining the macros
 * "CONTAINER_ID" and "N_ITEMS".
 *
 * Example:
 *  #define CONTAINER_ID NB_CONTAINER_QUEUE
 *  #define N_ITEMS 1000
 *  #include "iterator_TEMPLATE_UT.c"
 */
#include <stdbool.h>
#include <stdlib.h>
#include <stdint.h>

#include "nb/container_bot/container.h"
#include "nb/container_bot/iterator.h"

#include "test_library.h"
#include "test_add.h"

static bool check_create(void);
static bool check_set_container(void);
static bool check_clone(void);
static bool check_destroy(void);
static bool check_restart(void);
static bool check_get_next(void);
static bool check_has_more(void);
static bool check_has_more_with_1_item(void);

static vcn_container_t* get_container(int N);
static void insert_N_int32(vcn_container_t *cnt, int N);
static uint32_t keygen(const void *const val);
static bool are_equal(const void *const a, const void *const b);

inline int TEMPLATE_get_driver_id(void)
{
	return NB_DRIVER_UNIT_TEST;
}

void TEMPLATE_load_tests(void *tests_ptr)
{
	vcn_test_add(tests_ptr, check_create,
		     "Check create()");
	vcn_test_add(tests_ptr, check_set_container,
		     "Check set_container()");
	vcn_test_add(tests_ptr, check_clone,
		     "Check clone()");
	vcn_test_add(tests_ptr, check_destroy,
		     "Check destroy()");
	vcn_test_add(tests_ptr, check_restart,
		     "Check restart()");
	vcn_test_add(tests_ptr, check_get_next,
		     "Check get_next()");
	vcn_test_add(tests_ptr, check_has_more,
		     "Check has_more()");
	vcn_test_add(tests_ptr, check_has_more_with_1_item,
		     "Check has_more() with one item");
}

static bool check_create(void)
{
	vcn_iterator_t *iter = vcn_iterator_create();
	bool is_ok = (NULL != iter);
	if (is_ok)
		vcn_iterator_destroy(iter);
	return is_ok;
}

static bool check_set_container(void)
{
	vcn_container_t *cnt = get_container(N_ITEMS);
	vcn_iterator_t *iter = vcn_iterator_create();
	vcn_iterator_set_container(iter, cnt);
	bool is_ok = vcn_iterator_has_more(iter);
	vcn_iterator_destroy(iter);
	vcn_container_destroy(cnt);
	return is_ok;
}

static bool check_clone(void)
{
	vcn_container_t *cnt = get_container(N_ITEMS);
	vcn_iterator_t *iter = vcn_iterator_create();
	vcn_iterator_set_container(iter, cnt);
	vcn_iterator_t *cloned = vcn_iterator_clone(iter);
	bool is_ok = vcn_iterator_has_more(cloned);
	vcn_iterator_destroy(cloned);
	vcn_iterator_destroy(iter);
	vcn_container_destroy(cnt);
	return is_ok;	
}

static bool check_destroy(void)
{
	vcn_iterator_t *iter = vcn_iterator_create();
	vcn_iterator_destroy(iter);
	return true;
}

static bool check_restart(void)
{
	vcn_container_t *cnt = get_container(N_ITEMS);
	vcn_iterator_t *iter = vcn_iterator_create();
	vcn_iterator_set_container(iter, cnt);
	while (vcn_iterator_has_more(iter))
		vcn_iterator_get_next(iter);
	vcn_iterator_restart(iter);
	bool is_ok = vcn_iterator_has_more(iter);
	vcn_iterator_destroy(iter);
	vcn_container_destroy(cnt);
	return is_ok;
}

static bool check_get_next(void)
{
	vcn_container_t *cnt = get_container(N_ITEMS);
	vcn_iterator_t *iter = vcn_iterator_create();
	vcn_iterator_set_container(iter, cnt);
	bool is_ok = true;
	while (vcn_iterator_has_more(iter)) {
		const int32_t *val = vcn_iterator_get_next(iter);
		is_ok = is_ok && (NULL != val);
	}
	vcn_iterator_destroy(iter);
	vcn_container_destroy(cnt);
	return is_ok;
}

static bool check_has_more(void)
{
	vcn_container_t *cnt = get_container(N_ITEMS);
	vcn_iterator_t *iter = vcn_iterator_create();
	vcn_iterator_set_container(iter, cnt);
	int32_t counter = 0;
	while (vcn_iterator_has_more(iter)) {
		vcn_iterator_get_next(iter);
		counter += 1;
	}
	vcn_iterator_destroy(iter);
	vcn_container_destroy(cnt);
	return (N_ITEMS == counter);
}

static bool check_has_more_with_1_item(void)
{
	vcn_container_t *cnt = get_container(1);
	vcn_iterator_t *iter = vcn_iterator_create();
	vcn_iterator_set_container(iter, cnt);
	int32_t counter = 0;
	while (vcn_iterator_has_more(iter)) {
		vcn_iterator_get_next(iter);
		counter += 1;
	}
	vcn_iterator_destroy(iter);
	vcn_container_destroy(cnt);
	return (1 == counter);
}

static inline vcn_container_t* get_container(int N)
{
	vcn_container_t *cnt = vcn_container_create(CONTAINER_ID);
      	vcn_container_set_key_generator(cnt, keygen);
      	vcn_container_set_destroyer(cnt, free);
	insert_N_int32(cnt, N);
	return cnt;
}

static void insert_N_int32(vcn_container_t *cnt, int N)
{
	for (int32_t i = 0; i < N; i++) {
		int32_t *val = malloc(sizeof(*val));
		*val = i;
		vcn_container_insert(cnt, val);
	}
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
