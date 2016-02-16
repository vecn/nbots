/*
 * This file should not be compiled alone.
 * It must be included after defining the macro
 * "CONTAINER_ID".
 *
 * Example:
 *  #define CONTAINER_ID NB_QUEUE
 *  #include "container_TEMPLATE_CT.c"
 */

#include <stdbool.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "nb/container_bot/array.h"
#include "nb/container_bot/container.h"
#include "nb/statistics_bot.h"

#include "test_library.h"
#include "test_add.h"

static bool run_clone(void *data);
static bool run_merge(void *data);
static bool run_cast_to_QUEUE(void *data);
static bool run_cast_to_STACK(void *data);
static bool run_cast_to_SORTED(void *data);
static bool run_cast_to_HEAP(void *data);
static bool run_cast_to_HASH(void *data);
static bool run_clear(void *data);
static bool run_insert(void *data);
static bool run_get_first(void *data);
static bool run_delete_first(void *data);
static bool run_exist(void *data);
static bool run_delete(void *data);
static bool run_get_length(void *data);

static bool run_cast_to(void *data, int id);
static void* init_container_by_id(int id, void *init_data);
static uint32_t keygen(const void *const val);
static void* init_container(void *init_data);
static void destroy_container(void *data);
static void* clone_int8(const void *const a);
static void* init_cnt_and_clone(void *init_data);
static void destroy_cnt_and_clone(void *data);
static void* init_merge_with(void *init_data, int id);
static void destroy_merge_data(void *data);
static void* init_merge_with_QUEUE(void *init_data);
static void* init_merge_with_STACK(void *init_data);
static void* init_merge_with_SORTED(void *init_data);
static void* init_merge_with_HEAP(void *init_data);
static void* init_merge_with_HASH(void *init_data);
static void* init_cast_data(void *init_data);
static void destroy_cast_data(void *data);
static bool are_equal_int32(const void *const a, const void *const b);

inline int TEMPLATE_get_driver_id(void)
{
	return NB_DRIVER_COMPLEXITY_TEST;
}

void TEMPLATE_load_tests(void *tests_ptr)
{
	vcn_test_add_decomposed(tests_ptr, init_cnt_and_clone,
				run_clone, destroy_cnt_and_clone,
				"clone()");
	vcn_test_add_decomposed(tests_ptr, init_merge_with_QUEUE,
				run_merge, destroy_merge_data,
				"merge() with QUEUE");
	vcn_test_add_decomposed(tests_ptr, init_merge_with_STACK,
				run_merge, destroy_merge_data,
				"merge() with STACK");
	vcn_test_add_decomposed(tests_ptr, init_merge_with_SORTED,
				run_merge, destroy_merge_data,
				"merge() with SORTED");
	vcn_test_add_decomposed(tests_ptr, init_merge_with_HEAP,
				run_merge, destroy_merge_data,
				"merge() with HEAP");
	vcn_test_add_decomposed(tests_ptr, init_merge_with_HASH,
				run_merge, destroy_merge_data,
				"merge() with HASH");
	vcn_test_add_decomposed(tests_ptr, init_cast_data,
				run_cast_to_QUEUE, destroy_cast_data,
				"cast() to QUEUE");
	vcn_test_add_decomposed(tests_ptr, init_cast_data,
				run_cast_to_STACK, destroy_cast_data,
				"cast() to STACK");
	vcn_test_add_decomposed(tests_ptr, init_cast_data,
				run_cast_to_SORTED, destroy_cast_data,
				"cast() to SORTED");
	vcn_test_add_decomposed(tests_ptr, init_cast_data,
				run_cast_to_HEAP, destroy_cast_data,
				"cast() to HEAP");
	vcn_test_add_decomposed(tests_ptr, init_cast_data,
				run_cast_to_HASH, destroy_cast_data,
				"cast() to HASH");
	vcn_test_add_decomposed(tests_ptr, init_container,
				run_clear, destroy_container,
				"clear()");
	vcn_test_add_decomposed(tests_ptr, init_container,
				run_insert, destroy_container,
				"insert()");
	vcn_test_add_decomposed(tests_ptr, init_container,
				run_get_first, destroy_container,
				"get_first()");
	vcn_test_add_decomposed(tests_ptr, init_container,
				run_delete_first, destroy_container,
				"delete_first()");
	vcn_test_add_decomposed(tests_ptr, init_container,
				run_exist, destroy_container,
				"exist()");
	vcn_test_add_decomposed(tests_ptr, init_container,
				run_delete, destroy_container,
				"delete()");
	vcn_test_add_decomposed(tests_ptr, init_container,
				run_get_length, destroy_container,
				"get_length()");
}

static bool run_clone(void *data)
{
	void **containers = data;
	nb_container_set_cloner(containers[0], clone_int8);
	containers[1] = nb_container_clone(containers[0]);
	return true;
}

static bool run_merge(void *data)
{
	void **containers = data;
	nb_container_merge(containers[0], containers[1]);
	return true;
}

static inline bool run_cast_to_QUEUE(void *data)
{
	return run_cast_to(data, NB_QUEUE);
}

static inline bool run_cast_to_STACK(void *data)
{
	return run_cast_to(data, NB_STACK);
}

static inline bool run_cast_to_SORTED(void *data)
{
	return run_cast_to(data, NB_SORTED);
}

static inline bool run_cast_to_HEAP(void *data)
{
	return run_cast_to(data, NB_HEAP);
}

static inline bool run_cast_to_HASH(void *data)
{
	return run_cast_to(data, NB_HASH);
}

static bool run_clear(void *data)
{
	nb_container_clear(data);
	return true;
}

static bool run_insert(void *data)
{
	int32_t *val = malloc(sizeof(*val));
	*val = 1;
	nb_container_insert(data, val);
	return true;
}

static bool run_get_first(void *data)
{
	nb_container_get_first(data);
	return true;
}

static bool run_delete_first(void *data)
{
	void *val = nb_container_delete_first(data);
	bool is_ok = false;
	if (NULL != val) {
		free(val);
		is_ok = true;
	}
	return is_ok;
}

static bool run_exist(void *data)
{
	int32_t to_find = 1;
	nb_container_set_comparer(data, are_equal_int32);
	nb_container_exist(data, &to_find);
	return true;
}

static bool run_delete(void *data)
{
	int32_t to_delete = 1;
	nb_container_set_comparer(data, are_equal_int32);
	int32_t *val = nb_container_delete(data, &to_delete);
	bool is_ok = false;
	if (NULL != val) {
		if (to_delete == *val)
			is_ok = true;
		free(val);
	}
	return is_ok;
}

static bool run_get_length(void *data)
{
	nb_container_get_length(data);
	return true;
}

static bool run_cast_to(void *data, int id)
{
	nb_container_cast(data, id);
	return true;
}

static void* init_container_by_id(int id, void *init_data)
{
	int N = *((int*)init_data);
	nb_container_t *cnt = nb_container_create(id);
	nb_container_set_key_generator(cnt, keygen);
	nb_container_set_destroyer(cnt, free);
	int32_t seed = 1;
	for (int i = 0; i < N; i++) {
		int32_t *val = malloc(sizeof(*val));
		*val = seed;
		seed = vcn_statistics_lcg(seed);
		nb_container_insert(cnt, val);
	}
	return cnt;
}

static inline uint32_t keygen(const void *const val)
{
	return (uint32_t) *((int32_t*)val);
}

static inline void* init_container(void *init_data)
{
	return init_container_by_id(CONTAINER_ID, init_data);
}

static inline void destroy_container(void *data)
{
	nb_container_destroy(data);
}

static inline void* clone_int8(const void *const a)
{
	int8_t *b = malloc(sizeof(*b));
	*b = *((int8_t*)a);
	return b;
}

static void* init_cnt_and_clone(void *init_data)
{
	void **containers = calloc(2, sizeof(*containers));	
	containers[0] = init_container(init_data);
	return containers;
}

static void destroy_cnt_and_clone(void *data)
{
	void **containers = data;
	nb_container_destroy(containers[0]);
	nb_container_destroy(containers[1]);
	free(containers);
}

static void* init_merge_with(void *init_data, int id)
{
	void **containers = calloc(2, sizeof(*containers));
	containers[0] = init_container(init_data);
	containers[1] = init_container_by_id(id, init_data);
	return containers;
}

static void destroy_merge_data(void *data)
{
	void **containers = data;
	nb_container_set_destroyer(containers[0], free);
	nb_container_destroy(containers[0]);
	nb_container_destroy(containers[1]);
	free(containers);
}

static inline void* init_merge_with_QUEUE(void *init_data)
{
	return init_merge_with(init_data, NB_QUEUE);
}

static inline void* init_merge_with_STACK(void *init_data)
{
	return init_merge_with(init_data, NB_STACK);
}

static inline void* init_merge_with_SORTED(void *init_data)
{
	return init_merge_with(init_data, NB_SORTED);
}

static inline void* init_merge_with_HEAP(void *init_data)
{
	return init_merge_with(init_data, NB_HEAP);
}

static inline void* init_merge_with_HASH(void *init_data)
{
	return init_merge_with(init_data, NB_HASH);
}

static void* init_cast_data(void *init_data)
{
	return init_container(init_data);
}

static void destroy_cast_data(void *data)
{
	nb_container_destroy(data);
	free(data);
}

static inline bool are_equal_int32(const void *const a, const void *const b)
{
	return (0 == vcn_compare_int32(a, b));
}
