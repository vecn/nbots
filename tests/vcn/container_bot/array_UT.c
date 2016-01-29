#include <stdbool.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "vcn/container_bot/array.h"

#include "test_library.h"
#include "test_add.h"

#define PI (3.14159265358979323846)

typedef struct {
	float (*f)(float);
} data_t;

static bool check_array_get(void);
static bool check_swap_with_int8_ieqj(void);
static bool check_swap_with_int8_igtj(void);
static bool check_swap_with_int8_iltj(void);
static bool check_swap_with_int16(void);
static bool check_swap_with_int32(void);
static bool check_swap_with_int64(void);
static bool check_swap_with_uint8(void);
static bool check_swap_with_uint16(void);
static bool check_swap_with_uint32(void);
static bool check_swap_with_uint64(void);
static bool check_swap_with_float(void);
static bool check_swap_with_double(void);
static bool check_qsort(void);
static bool check_qsort_with_data(void);
static bool check_array_get_min(void);
static bool check_array_get_min_with_data(void);
static bool check_array_get_max(void);
static bool check_array_get_max_with_data(void);
static bool check_array_get_min_max(void);
static bool check_array_get_min_max_with_data(void);
static bool check_compare_char(void);
static bool check_compare_float(void);
static bool check_compare_double(void);
static bool check_compare_int8(void);
static bool check_compare_int16(void);
static bool check_compare_int32(void);
static bool check_compare_int64(void);
static bool check_compare_uint8(void);
static bool check_compare_uint16(void);
static bool check_compare_uint32(void);
static bool check_compare_uint64(void);

static int8_t compare_fnc(const void *const a, 
			  const void *const b,
			  const void *const data);

inline int vcn_test_get_driver_id(void)
{
	return VCN_DRIVER_UNIT_TEST;
}

void vcn_test_load_tests(void *tests_ptr)
{
	vcn_test_add(tests_ptr, check_array_get,
		     "Check array_get()");
	vcn_test_add(tests_ptr, check_swap_with_int8_ieqj,
		     "Check swap() with int8_t using i = j");
	vcn_test_add(tests_ptr, check_swap_with_int8_igtj,
		     "Check swap() with int8_t using i > j");
	vcn_test_add(tests_ptr, check_swap_with_int8_iltj,
		     "Check swap() with int8_t using j < i");
	vcn_test_add(tests_ptr, check_swap_with_int16,
		     "Check swap() with int16_t");
	vcn_test_add(tests_ptr, check_swap_with_int32,
		     "Check swap() with int32_t");
	vcn_test_add(tests_ptr, check_swap_with_int64,
		     "Check swap() with int64_t");
	vcn_test_add(tests_ptr, check_swap_with_uint8,
		     "Check swap() with uint8_t");
	vcn_test_add(tests_ptr, check_swap_with_uint16,
		     "Check swap() with uint16_t");
	vcn_test_add(tests_ptr, check_swap_with_uint32,
		     "Check swap() with uint32_t");
	vcn_test_add(tests_ptr, check_swap_with_uint64,
		     "Check swap() with uint64_t");
	vcn_test_add(tests_ptr, check_swap_with_float,
		     "Check swap() with float");
	vcn_test_add(tests_ptr, check_swap_with_double,
		     "Check swap() with double");
	vcn_test_add(tests_ptr, check_qsort,
		     "Check qsort()");
	vcn_test_add(tests_ptr, check_qsort_with_data,
		     "Check qsort() with extra data to compare");
	vcn_test_add(tests_ptr, check_array_get_min,
		     "Check array_get_min_id()");
	vcn_test_add(tests_ptr, check_array_get_min_with_data,
		     "Check array_get_min_id_wd()");
	vcn_test_add(tests_ptr, check_array_get_max,
		     "Check array_get_max_id()");
	vcn_test_add(tests_ptr, check_array_get_max_with_data,
		     "Check array_get_max_id_wd()");
	vcn_test_add(tests_ptr, check_array_get_min_max,
		     "Check array_get_min_max_ids()");
	vcn_test_add(tests_ptr, check_array_get_min_max_with_data,
		     "Check array_get_min_max_ids_wd()");
	vcn_test_add(tests_ptr, check_compare_char,
		     "Check compare_char()");
	vcn_test_add(tests_ptr, check_compare_float,
		     "Check compare_float()");
	vcn_test_add(tests_ptr, check_compare_double,
		     "Check compare_double()");
	vcn_test_add(tests_ptr, check_compare_int8,
		     "Check compare_int8()");
	vcn_test_add(tests_ptr, check_compare_int16,
		     "Check compare_int16()");
	vcn_test_add(tests_ptr, check_compare_int32,
		     "Check compare_int32()");
	vcn_test_add(tests_ptr, check_compare_int64,
		     "Check compare_int64()");
	vcn_test_add(tests_ptr, check_compare_uint8,
		     "Check compare_uint8()");
	vcn_test_add(tests_ptr, check_compare_uint16,
		     "Check compare_uint16()");
	vcn_test_add(tests_ptr, check_compare_uint32,
		     "Check compare_uint32()");
	vcn_test_add(tests_ptr, check_compare_uint64,
		     "Check compare_uint64()");
}

static bool check_array_get(void)
{
	uint8_t id = 2;
	int8_t array[5] = {2, 4, 6, 8, 10};
	int8_t *val = vcn_array_get(array, sizeof(*array), id);
	return  (array[id] == *val);
}

static bool check_swap_with_int8_ieqj(void)
{
	int8_t array[3] = {1, 2, 3};
	vcn_swap(array, 1, 1, sizeof(*array));
	return  (array[1] == 2);
}

static bool check_swap_with_int8_igtj(void)
{
	int8_t array[3] = {1, 2, 3};
	vcn_swap(array, 2, 0, sizeof(*array));
	return  (array[0] == 3) && (array[2] == 1);
}

static bool check_swap_with_int8_iltj(void)
{
	int8_t array[3] = {1, 2, 3};
	vcn_swap(array, 0, 2, sizeof(*array));
	return  (array[0] == 3) && (array[2] == 1);
}

static bool check_swap_with_int16(void)
{
	int16_t array[3] = {1, 2, 3};
	vcn_swap(array, 0, 2, sizeof(*array));
	return  (array[0] == 3) && (array[2] == 1);
}

static bool check_swap_with_int32(void)
{
	int32_t array[3] = {1, 2, 3};
	vcn_swap(array, 0, 2, sizeof(*array));
	return  (array[0] == 3) && (array[2] == 1);
}

static bool check_swap_with_int64(void)
{
	int64_t array[3] = {1, 2, 3};
	vcn_swap(array, 0, 2, sizeof(*array));
	return  (array[0] == 3) && (array[2] == 1);
}

static bool check_swap_with_uint8(void)
{
	uint8_t array[3] = {1U, 2U, 3U};
	vcn_swap(array, 0, 2, sizeof(*array));
	return  (array[0] == 3) && (array[2] == 1);
}

static bool check_swap_with_uint16(void)
{
	uint16_t array[3] = {1U, 2U, 3U};
	vcn_swap(array, 0, 2, sizeof(*array));
	return  (array[0] == 3) && (array[2] == 1);
}

static bool check_swap_with_uint32(void)
{
	uint32_t array[3] = {1U, 2U, 3U};
	vcn_swap(array, 0, 2, sizeof(*array));
	return  (array[0] == 3) && (array[2] == 1);
}

static bool check_swap_with_uint64(void)
{
	uint64_t array[3] = {1U, 2U, 3U};
	vcn_swap(array, 0, 2, sizeof(*array));
	return  (array[0] == 3) && (array[2] == 1);
}

static bool check_swap_with_float(void)
{
	float array[3] = {1.0f, 2.0f, 3.0f};
	vcn_swap(array, 0, 2, sizeof(*array));
	return  (array[0] == 3) && (array[2] == 1);
}

static bool check_swap_with_double(void)
{
	double array[3] = {1.0, 2.0, 3.0};
	vcn_swap(array, 0, 2, sizeof(*array));
	return  (array[0] == 3) && (array[2] == 1);
}

static bool check_qsort(void)
{
	uint32_t N = 100;
	int8_t *array = calloc(N, sizeof(*array));
	for (int i = 0; i < N; i++)
		array[i] = N - i;
	vcn_qsort(array, N, sizeof(*array), vcn_compare_int8);
	bool are_sorted = true;
	for (int i = 1; i < N; i++) {
		if (array[i-1] > array[i]) {
			are_sorted = false;
			break;
		}
	}
	free(array);
	return are_sorted;
}

static bool check_qsort_with_data(void)
{
	uint32_t N = 100;
	float *array = calloc(N, sizeof(*array));
	for (int i = 0; i < N; i++)
		array[i] = 2.0f * PI * (float)i/(N - 1.0f);
	data_t data;
	data.f = sinf;

	vcn_qsort_wd(array, N, sizeof(*array), compare_fnc, &data);
	bool are_sorted = true;
	for (int i = 1; i < N; i++) {
		if (data.f(array[i-1]) > data.f(array[i])) {
			are_sorted = false;
			break;
		}
	}
	free(array);
	return are_sorted;
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
	uint32_t N = 100;
	int8_t *array = calloc(N, sizeof(*array));
	for (int i = 0; i < N; i++)
		array[i] = N - i;

	uint32_t imin = vcn_array_get_min_id(array, N, sizeof(*array),
					     vcn_compare_int8);
	bool success = (1 == array[imin]);
	free(array);
	return success;
}

static bool check_array_get_min_with_data(void)
{
	uint32_t N = 100;
	float *array = calloc(N, sizeof(*array));
	for (int i = 0; i < N; i++)
		array[i] = (float)(N - i);
	data_t data;
	data.f = sqrtf;

	uint32_t imin = vcn_array_get_min_id_wd(array, N, sizeof(*array),
						compare_fnc, &data);
	bool success = fabsf(data.f(1.0f) - data.f(array[imin])) < 1e-5;
	free(array);
	return success;
}

static bool check_array_get_max(void)
{
	uint32_t N = 100;
	int8_t *array = calloc(N, sizeof(*array));
	for (int i = 0; i < N; i++)
		array[i] = N - i;

	uint32_t imax = vcn_array_get_max_id(array, N, sizeof(*array),
					     vcn_compare_int8);
	bool success = (N == array[imax]);
	free(array);
	return success;
}

static bool check_array_get_max_with_data(void)
{
	uint32_t N = 100;
	float *array = calloc(N, sizeof(*array));
	for (int i = 0; i < N; i++)
		array[i] = (float)(N - i);
	data_t data;
	data.f = sqrtf;

	uint32_t imax = vcn_array_get_max_id_wd(array, N, sizeof(*array),
						compare_fnc, &data);
	bool success = fabsf(data.f((float)N) - data.f(array[imax])) < 1e-5;
	free(array);
	return success;
}

static bool check_array_get_min_max(void)
{
	uint32_t N = 100;
	int8_t *array = calloc(N, sizeof(*array));
	for (int i = 0; i < N; i++)
		array[i] = N - i;

	uint32_t imin;
	uint32_t imax;
	vcn_array_get_min_max_ids(array, N, sizeof(*array), vcn_compare_int8,
				  &imin, &imax);
	bool success = (1 == array[imin]) && (N == array[imax]);
	free(array);
	return success;
}

static bool check_array_get_min_max_with_data(void)
{
	uint32_t N = 100;
	float *array = calloc(N, sizeof(*array));
	for (int i = 0; i < N; i++)
		array[i] = (float)(N - i);
	data_t data;
	data.f = sqrtf;

	uint32_t imin;
	uint32_t imax;
	vcn_array_get_min_max_ids_wd(array, N, sizeof(*array),
				     compare_fnc, &data, &imin, &imax);
	bool success = 
		(fabsf(data.f((float)1) - data.f(array[imin])) < 1e-5) &&
		(fabsf(data.f((float)N) - data.f(array[imax])) < 1e-5);
	free(array);
	return success;
}

static bool check_compare_char(void)
{
	char a = 0;
	char b = 0;
	char c = 1;
	int8_t cmp_eq = vcn_compare_char(&a, &b);
	int8_t cmp_lt = vcn_compare_char(&b, &c);
	int8_t cmp_gt = vcn_compare_char(&c, &a);
	return (cmp_eq == 0) && (cmp_lt == -1) && (cmp_gt == 1);
}

static bool check_compare_float(void)
{
	float a = 0.0f;
	float b = 0.0f;
	float c = 1.0f;
	int8_t cmp_eq = vcn_compare_float(&a, &b);
	int8_t cmp_lt = vcn_compare_float(&b, &c);
	int8_t cmp_gt = vcn_compare_float(&c, &a);
	return (cmp_eq == 0) && (cmp_lt == -1) && (cmp_gt == 1);
}

static bool check_compare_double(void)
{
	double a = 0.0;
	double b = 0.0;
	double c = 1.0;
	int8_t cmp_eq = vcn_compare_double(&a, &b);
	int8_t cmp_lt = vcn_compare_double(&b, &c);
	int8_t cmp_gt = vcn_compare_double(&c, &a);
	return (cmp_eq == 0) && (cmp_lt == -1) && (cmp_gt == 1);
}

static bool check_compare_int8(void)
{
	int8_t a = 0;
	int8_t b = 0;
	int8_t c = 1;
	int8_t cmp_eq = vcn_compare_int8(&a, &b);
	int8_t cmp_lt = vcn_compare_int8(&b, &c);
	int8_t cmp_gt = vcn_compare_int8(&c, &a);
	return (cmp_eq == 0) && (cmp_lt == -1) && (cmp_gt == 1);
}

static bool check_compare_int16(void)
{
	int16_t a = 0;
	int16_t b = 0;
	int16_t c = 1;
	int8_t cmp_eq = vcn_compare_int16(&a, &b);
	int8_t cmp_lt = vcn_compare_int16(&b, &c);
	int8_t cmp_gt = vcn_compare_int16(&c, &a);
	return (cmp_eq == 0) && (cmp_lt == -1) && (cmp_gt == 1);
}

static bool check_compare_int32(void)
{
	int32_t a = 0;
	int32_t b = 0;
	int32_t c = 1;
	int8_t cmp_eq = vcn_compare_int32(&a, &b);
	int8_t cmp_lt = vcn_compare_int32(&b, &c);
	int8_t cmp_gt = vcn_compare_int32(&c, &a);
	return (cmp_eq == 0) && (cmp_lt == -1) && (cmp_gt == 1);
}

static bool check_compare_int64(void)
{
	int64_t a = 0;
	int64_t b = 0;
	int64_t c = 1;
	int8_t cmp_eq = vcn_compare_int64(&a, &b);
	int8_t cmp_lt = vcn_compare_int64(&b, &c);
	int8_t cmp_gt = vcn_compare_int64(&c, &a);
	return (cmp_eq == 0) && (cmp_lt == -1) && (cmp_gt == 1);
}

static bool check_compare_uint8(void)
{
	uint8_t a = 0U;
	uint8_t b = 0U;
	uint8_t c = 1U;
	int8_t cmp_eq = vcn_compare_uint8(&a, &b);
	int8_t cmp_lt = vcn_compare_uint8(&b, &c);
	int8_t cmp_gt = vcn_compare_uint8(&c, &a);
	return (cmp_eq == 0) && (cmp_lt == -1) && (cmp_gt == 1);
}

static bool check_compare_uint16(void)
{
	uint16_t a = 0U;
	uint16_t b = 0U;
	uint16_t c = 1U;
	int8_t cmp_eq = vcn_compare_uint16(&a, &b);
	int8_t cmp_lt = vcn_compare_uint16(&b, &c);
	int8_t cmp_gt = vcn_compare_uint16(&c, &a);
	return (cmp_eq == 0) && (cmp_lt == -1) && (cmp_gt == 1);
}

static bool check_compare_uint32(void)
{
	uint32_t a = 0U;
	uint32_t b = 0U;
	uint32_t c = 1U;
	int8_t cmp_eq = vcn_compare_uint32(&a, &b);
	int8_t cmp_lt = vcn_compare_uint32(&b, &c);
	int8_t cmp_gt = vcn_compare_uint32(&c, &a);
	return (cmp_eq == 0) && (cmp_lt == -1) && (cmp_gt == 1);
}

static bool check_compare_uint64(void)
{
	uint64_t a = 0U;
	uint64_t b = 0U;
	uint64_t c = 1U;
	int8_t cmp_eq = vcn_compare_uint64(&a, &b);
	int8_t cmp_lt = vcn_compare_uint64(&b, &c);
	int8_t cmp_gt = vcn_compare_uint64(&c, &a);
	return (cmp_eq == 0) && (cmp_lt == -1) && (cmp_gt == 1);
}
