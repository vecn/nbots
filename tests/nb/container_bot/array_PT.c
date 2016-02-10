#include <stdbool.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "nb/container_bot/array.h"

#include "test_library.h"
#include "test_add.h"

#define ARRAY_SIZE 100000

typedef struct {
	float (*f)(float);
} data_t;

static bool check_swap(void);
static bool check_qsort(void);
static bool check_qsort_wd(void);
static bool check_array_get_min(void);
static bool check_array_get_min_wd(void);
static bool check_array_get_max(void);
static bool check_array_get_max_wd(void);
static bool check_array_get_min_max(void);
static bool check_array_get_min_max_wd(void);

static void init_int8_array(uint32_t N, int8_t array[]);
static void init_float_array(uint32_t N, float array[]);
static int8_t compare_fnc(const void *const a, 
			  const void *const b,
			  const void *const data);

inline int vcn_test_get_driver_id(void)
{
	return VCN_DRIVER_PERFORMANCE_TEST;
}

void vcn_test_load_tests(void *tests_ptr)
{
	vcn_test_add(tests_ptr, check_swap,
		     "swap()");
	vcn_test_add(tests_ptr, check_qsort,
		     "qsort()");
	vcn_test_add(tests_ptr, check_qsort_wd,
		     "qsort_wd()");
	vcn_test_add(tests_ptr, check_array_get_min,
		     "array_get_min_id()");
	vcn_test_add(tests_ptr, check_array_get_min_wd,
		     "array_get_min_id_wd()");
	vcn_test_add(tests_ptr, check_array_get_max,
		     "array_get_max_id()");
	vcn_test_add(tests_ptr, check_array_get_max_wd,
		     "array_get_max_id_wd()");
	vcn_test_add(tests_ptr, check_array_get_min_max,
		     "array_get_min_max_ids()");
	vcn_test_add(tests_ptr, check_array_get_min_max_wd,
		     "array_get_min_max_ids_wd()");
}

static bool check_swap(void)
{
	char array[2] = {0, 1};
	vcn_swap(array, 0, 1, sizeof(*array));
	return true;
}

static bool check_qsort(void)
{
	int8_t array[ARRAY_SIZE];
	init_int8_array(ARRAY_SIZE, array);
	vcn_qsort(array, ARRAY_SIZE, sizeof(*array), vcn_compare_int8);
	return true;
}

static void init_int8_array(uint32_t N, int8_t array[])
{
	for (uint32_t i = 0; i < N; i++)
		array[i] = N - i;
}

static bool check_qsort_wd(void)
{
	float array[ARRAY_SIZE];
	init_float_array(ARRAY_SIZE, array);
	data_t data;
	data.f = sinf;
	vcn_qsort_wd(array, ARRAY_SIZE, sizeof(*array), compare_fnc, &data);
	return true;
}

static void init_float_array(uint32_t N, float array[])
{
	for (uint32_t i = 0; i < N; i++)
		array[i] = (float) (N - i);
}

static int8_t compare_fnc(const void *const restrict a, 
			  const void *const restrict b,
			  const void *const data)
{
	float a_val = *((float*)a);
	float b_val = *((float*)b);
	const data_t *const dat = data;
	a_val = dat->f(a_val);
	b_val = dat->f(b_val);
	int8_t out = 0;
	if (a_val > b_val)
		out = 1;
	else if (a_val < b_val)
		out = -1;
	return out;
}

static bool check_array_get_min(void)
{
	int8_t array[ARRAY_SIZE];
	init_int8_array(ARRAY_SIZE, array);
	vcn_array_get_min_id(array, ARRAY_SIZE, sizeof(*array),
			     vcn_compare_int8);
	return true;
}

static bool check_array_get_min_wd(void)
{
	float array[ARRAY_SIZE];
	init_float_array(ARRAY_SIZE, array);
	data_t data;
	data.f = sinf;
	vcn_array_get_min_id_wd(array, ARRAY_SIZE, sizeof(*array), compare_fnc,
				&data);
	return true;
}

static bool check_array_get_max(void)
{
	int8_t array[ARRAY_SIZE];
	init_int8_array(ARRAY_SIZE, array);
	vcn_array_get_max_id(array, ARRAY_SIZE, sizeof(*array),
			     vcn_compare_int8);
	return true;
}

static bool check_array_get_max_wd(void)
{
	float array[ARRAY_SIZE];
	init_float_array(ARRAY_SIZE, array);
	data_t data;
	data.f = sinf;
	vcn_array_get_max_id_wd(array, ARRAY_SIZE, sizeof(*array), compare_fnc,
				&data);
	return true;
}

static bool check_array_get_min_max(void)
{
	int8_t array[ARRAY_SIZE];
	init_int8_array(ARRAY_SIZE, array);
	uint32_t imin;
	uint32_t imax;
	vcn_array_get_min_max_ids(array, ARRAY_SIZE, sizeof(*array),
				  vcn_compare_int8, &imin, &imax);
	return true;
}

static bool check_array_get_min_max_wd(void)
{
	float array[ARRAY_SIZE];
	init_float_array(ARRAY_SIZE, array);
	data_t data;
	data.f = sinf;
	uint32_t imin;
	uint32_t imax;
	vcn_array_get_min_max_ids_wd(array, ARRAY_SIZE, sizeof(*array),
				     compare_fnc, &data, &imin, &imax);
	return true;
}
