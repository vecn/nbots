#include <stdbool.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "nb/container_bot/array.h"

#include "test_library.h"
#include "test_add.h"

typedef struct {
	int32_t N;
	int16_t type_size;
	void *base;
} array_t;

typedef struct {
	float (*f)(float);
} fnc_t;

static bool run_qsort(void *data);
static bool run_qsort_wd(void *data);
static bool run_array_get_min(void *data);
static bool run_array_get_min_wd(void *data);
static bool run_array_get_max(void *data);
static bool run_array_get_max_wd(void *data);
static bool run_array_get_min_max(void *data);
static bool run_array_get_min_max_wd(void *data);

static void* init_int8_array(void *init_data);
static void destroy_array(void *data);
static void* init_float_array(void *init_data);
static int8_t compare_fnc(const void *const a, 
			  const void *const b,
			  const void *const data);

inline int vcn_test_get_driver_id(void)
{
	return VCN_DRIVER_COMPLEXITY_TEST;
}

void vcn_test_load_tests(void *tests_ptr)
{
	vcn_test_add_decomposed(tests_ptr, init_int8_array,
				run_qsort, destroy_array,
				"qsort()");
	vcn_test_add_decomposed(tests_ptr, init_float_array,
				run_qsort_wd, destroy_array,
				"qsort_wd()");
	vcn_test_add_decomposed(tests_ptr, init_int8_array,
				run_array_get_min, destroy_array,
				"array_get_min_id()");
	vcn_test_add_decomposed(tests_ptr, init_float_array,
				run_array_get_min_wd, destroy_array,
				"array_get_min_id_wd()");
	vcn_test_add_decomposed(tests_ptr, init_int8_array,
				run_array_get_max, destroy_array,
				"array_get_max_id()");
	vcn_test_add_decomposed(tests_ptr, init_float_array,
				run_array_get_max_wd, destroy_array,
				"array_get_max_id_wd()");
	vcn_test_add_decomposed(tests_ptr, init_int8_array,
				run_array_get_min_max, destroy_array,
				"array_get_min_max_ids()");
	vcn_test_add_decomposed(tests_ptr, init_float_array,
				run_array_get_min_max_wd, destroy_array,
				"array_get_min_max_ids_wd()");
}


static bool run_qsort(void *data)
{
	array_t *array = data;
	vcn_qsort(array->base, array->N, 
		  array->type_size, vcn_compare_int8);
	return true;
}

static bool run_qsort_wd(void *data)
{
	array_t *array = data;
	fnc_t fnc;
	fnc.f = sinf;
	vcn_qsort_wd(array->base, array->N, 
		     array->type_size, compare_fnc, &fnc);
	return true;
}

static bool run_array_get_min(void *data)
{
	array_t *array = data;
	vcn_array_get_min_id(array->base, array->N, 
			     array->type_size, vcn_compare_int8);
	return true;
}

static bool run_array_get_min_wd(void *data)
{
	array_t *array = data;
	fnc_t fnc;
	fnc.f = sinf;
	vcn_array_get_min_id_wd(array->base, array->N, 
				array->type_size, compare_fnc, &fnc);
	return true;
}

static bool run_array_get_max(void *data)
{
	array_t *array = data;
	vcn_array_get_max_id(array->base, array->N, 
			     array->type_size, vcn_compare_int8);
	return true;
}

static bool run_array_get_max_wd(void *data)
{
	array_t *array = data;
	fnc_t fnc;
	fnc.f = sinf;
	vcn_array_get_max_id_wd(array->base, array->N, 
				array->type_size, compare_fnc, &fnc);
	return true;
}

static bool run_array_get_min_max(void *data)
{
	array_t *array = data;
	uint32_t imin;
	uint32_t imax;
	vcn_array_get_min_max_ids(array->base, array->N, 
				  array->type_size, vcn_compare_int8, 
				  &imin, &imax);
	return true;
}

static bool run_array_get_min_max_wd(void *data)
{
	array_t *array = data;
	uint32_t imin;
	uint32_t imax;
	fnc_t fnc;
	fnc.f = sinf;
	vcn_array_get_min_max_ids_wd(array->base, array->N, 
				     array->type_size, compare_fnc, &fnc,
				     &imin, &imax);
	return true;
}

static void* init_int8_array(void *init_data)
{
	
	array_t *array = malloc(sizeof(array_t));
	array->N = *((int*)init_data);
	array->type_size = sizeof(int8_t);
	int8_t *base = malloc(array->N * array->type_size);
	for (int i = 0; i < array->N; i++)
		base[i] = (int8_t)(100 * sinf((float)(array->N - i)));
	array->base = base;
	return array;
}

static void destroy_array(void *data)
{
	array_t *array = data;
	free(array->base);
	free(array);
}

static void* init_float_array(void *init_data)
{
	array_t *array = malloc(sizeof(array_t));
	array->N = *((int*)init_data);
	array->type_size = sizeof(float);
	float *base = malloc(array->N * array->type_size);
	for (int i = 0; i < array->N; i++)
		base[i] = cosf((float)(array->N - i));
	array->base = base;
	return array;
}

static int8_t compare_fnc(const void *const restrict a, 
			  const void *const restrict b,
			  const void *const data)
{
	float a_val = *((float*)a);
	float b_val = *((float*)b);
	const fnc_t *const fnc = data;
	a_val = fnc->f(a_val);
	b_val = fnc->f(b_val);
	return vcn_compare_float(&a_val, &b_val);
}
