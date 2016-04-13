#include <stdbool.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include <CUnit/Basic.h>

#include "nb/container_bot/array.h"

#define PI (3.14159265358979323846)

typedef struct {
	float (*f)(float);
} data_t;

static int suite_init(void);
static int suite_clean(void);

static void test_array_get(void);
static void test_swap_with_int8_ieqj(void);
static void test_swap_with_int8_igtj(void);
static void test_swap_with_int8_iltj(void);
static void test_swap_with_int16(void);
static void test_swap_with_int32(void);
static void test_swap_with_int64(void);
static void test_swap_with_uint8(void);
static void test_swap_with_uint16(void);
static void test_swap_with_uint32(void);
static void test_swap_with_uint64(void);
static void test_swap_with_float(void);
static void test_swap_with_double(void);
static void test_qsort(void);
static void test_qsort_with_data(void);
static void test_array_get_min(void);
static void test_array_get_min_with_data(void);
static void test_array_get_max(void);
static void test_array_get_max_with_data(void);
static void test_array_get_min_max(void);
static void test_array_get_min_max_with_data(void);
static void test_compare_char(void);
static void test_compare_float(void);
static void test_compare_double(void);
static void test_compare_int8(void);
static void test_compare_int16(void);
static void test_compare_int32(void);
static void test_compare_int64(void);
static void test_compare_uint8(void);
static void test_compare_uint16(void);
static void test_compare_uint32(void);
static void test_compare_uint64(void);

static int8_t compare_fnc(const void *const a, 
			  const void *const b,
			  const void *const data);

void cunit_nb_container_bot_array(void)
{
	CU_pSuite suite = CU_add_suite("nb/container_bot/array.c",
				       suite_init, suite_clean);

	CU_add_test(suite, "array_get()", test_array_get);
	CU_add_test(suite, "swap() with int8_t using i eq j",
		    test_swap_with_int8_ieqj);
	CU_add_test(suite, "swap() with int8_t using i gt j",
		    test_swap_with_int8_igtj);
	CU_add_test(suite, "swap() with int8_t using j lt i",
		    test_swap_with_int8_iltj);
	CU_add_test(suite, "swap() with int16_t",
		    test_swap_with_int16);
	CU_add_test(suite, "swap() with int32_t",
		    test_swap_with_int32);
	CU_add_test(suite, "swap() with int64_t",
		    test_swap_with_int64);
	CU_add_test(suite, "swap() with uint8_t",
		    test_swap_with_uint8);
	CU_add_test(suite, "swap() with uint16_t",
		    test_swap_with_uint16);
	CU_add_test(suite, "swap() with uint32_t",
		    test_swap_with_uint32);
	CU_add_test(suite, "swap() with uint64_t",
		    test_swap_with_uint64);
	CU_add_test(suite, "swap() with float",
		    test_swap_with_float);
	CU_add_test(suite, "swap() with double",
		    test_swap_with_double);
	CU_add_test(suite, "qsort()", test_qsort);
	CU_add_test(suite, "qsort() with extra data to compare",
		    test_qsort_with_data);
	CU_add_test(suite, "array_get_min_id()",
		    test_array_get_min);
	CU_add_test(suite, "array_get_min_id_wd()",
		    test_array_get_min_with_data);
	CU_add_test(suite, "array_get_max_id()",
		    test_array_get_max);
	CU_add_test(suite, "array_get_max_id_wd()",
		    test_array_get_max_with_data);
	CU_add_test(suite, "array_get_min_max_ids()",
		    test_array_get_min_max);
	CU_add_test(suite, "array_get_min_max_ids_wd()",
		    test_array_get_min_max_with_data);
	CU_add_test(suite, "compare_char()", test_compare_char);
	CU_add_test(suite, "compare_float()", test_compare_float);
	CU_add_test(suite, "compare_double()", test_compare_double);
	CU_add_test(suite, "compare_int8()", test_compare_int8);
	CU_add_test(suite, "compare_int16()", test_compare_int16);
	CU_add_test(suite, "compare_int32()", test_compare_int32);
	CU_add_test(suite, "compare_int64()", test_compare_int64);
	CU_add_test(suite, "compare_uint8()", test_compare_uint8);
	CU_add_test(suite, "compare_uint16()", test_compare_uint16);
	CU_add_test(suite, "compare_uint32()", test_compare_uint32);
	CU_add_test(suite, "compare_uint64()", test_compare_uint64);
}

static int suite_init(void)
{
	return 0;
}

static int suite_clean(void)
{
	return 0;
}

static void test_array_get(void)
{
	uint8_t id = 2;
	int8_t array[5] = {2, 4, 6, 8, 10};
	int8_t *val = vcn_array_get(array, sizeof(*array), id);
	CU_ASSERT(array[id] == *val);
}

static void test_swap_with_int8_ieqj(void)
{
	int8_t array[3] = {1, 2, 3};
	vcn_swap(array, 1, 1, sizeof(*array));
	CU_ASSERT(array[1] == 2);
}

static void test_swap_with_int8_igtj(void)
{
	int8_t array[3] = {1, 2, 3};
	vcn_swap(array, 2, 0, sizeof(*array));
	CU_ASSERT(array[0] == 3);
	CU_ASSERT(array[2] == 1);
}

static void test_swap_with_int8_iltj(void)
{
	int8_t array[3] = {1, 2, 3};
	vcn_swap(array, 0, 2, sizeof(*array));
	CU_ASSERT(array[0] == 3);
	CU_ASSERT(array[2] == 1);
}

static void test_swap_with_int16(void)
{
	int16_t array[3] = {1, 2, 3};
	vcn_swap(array, 0, 2, sizeof(*array));
	CU_ASSERT(array[0] == 3);
	CU_ASSERT(array[2] == 1);
}

static void test_swap_with_int32(void)
{
	int32_t array[3] = {1, 2, 3};
	vcn_swap(array, 0, 2, sizeof(*array));
	CU_ASSERT(array[0] == 3);
	CU_ASSERT(array[2] == 1);
}

static void test_swap_with_int64(void)
{
	int64_t array[3] = {1, 2, 3};
	vcn_swap(array, 0, 2, sizeof(*array));
	CU_ASSERT(array[0] == 3);
	CU_ASSERT(array[2] == 1);
}

static void test_swap_with_uint8(void)
{
	uint8_t array[3] = {1U, 2U, 3U};
	vcn_swap(array, 0, 2, sizeof(*array));
	CU_ASSERT(array[0] == 3);
	CU_ASSERT(array[2] == 1);
}

static void test_swap_with_uint16(void)
{
	uint16_t array[3] = {1U, 2U, 3U};
	vcn_swap(array, 0, 2, sizeof(*array));
	CU_ASSERT(array[0] == 3);
	CU_ASSERT(array[2] == 1);
}

static void test_swap_with_uint32(void)
{
	uint32_t array[3] = {1U, 2U, 3U};
	vcn_swap(array, 0, 2, sizeof(*array));
	CU_ASSERT(array[0] == 3);
	CU_ASSERT(array[2] == 1);
}

static void test_swap_with_uint64(void)
{
	uint64_t array[3] = {1U, 2U, 3U};
	vcn_swap(array, 0, 2, sizeof(*array));
	CU_ASSERT(array[0] == 3);
	CU_ASSERT(array[2] == 1);
}

static void test_swap_with_float(void)
{
	float array[3] = {1.0f, 2.0f, 3.0f};
	vcn_swap(array, 0, 2, sizeof(*array));
	CU_ASSERT(array[0] == 3);
	CU_ASSERT(array[2] == 1);
}

static void test_swap_with_double(void)
{
	double array[3] = {1.0, 2.0, 3.0};
	vcn_swap(array, 0, 2, sizeof(*array));
	CU_ASSERT(array[0] == 3);
	CU_ASSERT(array[2] == 1);
}

static void test_qsort(void)
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
	CU_ASSERT(are_sorted);
}

static void test_qsort_with_data(void)
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
	CU_ASSERT(are_sorted);
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

static void test_array_get_min(void)
{
	uint32_t N = 100;
	int8_t *array = calloc(N, sizeof(*array));
	for (int i = 0; i < N; i++)
		array[i] = N - i;

	uint32_t imin = vcn_array_get_min_id(array, N, sizeof(*array),
					     vcn_compare_int8);
	bool success = (1 == array[imin]);
	free(array);
	CU_ASSERT(success);
}

static void test_array_get_min_with_data(void)
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
	CU_ASSERT(success);
}

static void test_array_get_max(void)
{
	uint32_t N = 100;
	int8_t *array = calloc(N, sizeof(*array));
	for (int i = 0; i < N; i++)
		array[i] = N - i;

	uint32_t imax = vcn_array_get_max_id(array, N, sizeof(*array),
					     vcn_compare_int8);
	bool success = (N == array[imax]);
	free(array);
	CU_ASSERT(success);
}

static void test_array_get_max_with_data(void)
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
	CU_ASSERT(success);
}

static void test_array_get_min_max(void)
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
	CU_ASSERT(success);
}

static void test_array_get_min_max_with_data(void)
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
	CU_ASSERT(success);
}

static void test_compare_char(void)
{
	char a = 0;
	char b = 0;
	char c = 1;
	int8_t cmp_eq = vcn_compare_char(&a, &b);
	int8_t cmp_lt = vcn_compare_char(&b, &c);
	int8_t cmp_gt = vcn_compare_char(&c, &a);
	CU_ASSERT(cmp_eq == 0);
	CU_ASSERT(cmp_lt == -1);
	CU_ASSERT(cmp_gt == 1);
}

static void test_compare_float(void)
{
	float a = 0.0f;
	float b = 0.0f;
	float c = 1.0f;
	int8_t cmp_eq = vcn_compare_float(&a, &b);
	int8_t cmp_lt = vcn_compare_float(&b, &c);
	int8_t cmp_gt = vcn_compare_float(&c, &a);
	CU_ASSERT(cmp_eq == 0);
	CU_ASSERT(cmp_lt == -1);
	CU_ASSERT(cmp_gt == 1);
}

static void test_compare_double(void)
{
	double a = 0.0;
	double b = 0.0;
	double c = 1.0;
	int8_t cmp_eq = vcn_compare_double(&a, &b);
	int8_t cmp_lt = vcn_compare_double(&b, &c);
	int8_t cmp_gt = vcn_compare_double(&c, &a);
	CU_ASSERT(cmp_eq == 0);
	CU_ASSERT(cmp_lt == -1);
	CU_ASSERT(cmp_gt == 1);
}

static void test_compare_int8(void)
{
	int8_t a = 0;
	int8_t b = 0;
	int8_t c = 1;
	int8_t cmp_eq = vcn_compare_int8(&a, &b);
	int8_t cmp_lt = vcn_compare_int8(&b, &c);
	int8_t cmp_gt = vcn_compare_int8(&c, &a);
	CU_ASSERT(cmp_eq == 0);
	CU_ASSERT(cmp_lt == -1);
	CU_ASSERT(cmp_gt == 1);
}

static void test_compare_int16(void)
{
	int16_t a = 0;
	int16_t b = 0;
	int16_t c = 1;
	int8_t cmp_eq = vcn_compare_int16(&a, &b);
	int8_t cmp_lt = vcn_compare_int16(&b, &c);
	int8_t cmp_gt = vcn_compare_int16(&c, &a);
	CU_ASSERT(cmp_eq == 0);
	CU_ASSERT(cmp_lt == -1);
	CU_ASSERT(cmp_gt == 1);
}

static void test_compare_int32(void)
{
	int32_t a = 0;
	int32_t b = 0;
	int32_t c = 1;
	int8_t cmp_eq = vcn_compare_int32(&a, &b);
	int8_t cmp_lt = vcn_compare_int32(&b, &c);
	int8_t cmp_gt = vcn_compare_int32(&c, &a);
	CU_ASSERT(cmp_eq == 0);
	CU_ASSERT(cmp_lt == -1);
	CU_ASSERT(cmp_gt == 1);
}

static void test_compare_int64(void)
{
	int64_t a = 0;
	int64_t b = 0;
	int64_t c = 1;
	int8_t cmp_eq = vcn_compare_int64(&a, &b);
	int8_t cmp_lt = vcn_compare_int64(&b, &c);
	int8_t cmp_gt = vcn_compare_int64(&c, &a);
	CU_ASSERT(cmp_eq == 0);
	CU_ASSERT(cmp_lt == -1);
	CU_ASSERT(cmp_gt == 1);
}

static void test_compare_uint8(void)
{
	uint8_t a = 0U;
	uint8_t b = 0U;
	uint8_t c = 1U;
	int8_t cmp_eq = vcn_compare_uint8(&a, &b);
	int8_t cmp_lt = vcn_compare_uint8(&b, &c);
	int8_t cmp_gt = vcn_compare_uint8(&c, &a);
	CU_ASSERT(cmp_eq == 0);
	CU_ASSERT(cmp_lt == -1);
	CU_ASSERT(cmp_gt == 1);
}

static void test_compare_uint16(void)
{
	uint16_t a = 0U;
	uint16_t b = 0U;
	uint16_t c = 1U;
	int8_t cmp_eq = vcn_compare_uint16(&a, &b);
	int8_t cmp_lt = vcn_compare_uint16(&b, &c);
	int8_t cmp_gt = vcn_compare_uint16(&c, &a);
	CU_ASSERT(cmp_eq == 0);
	CU_ASSERT(cmp_lt == -1);
	CU_ASSERT(cmp_gt == 1);
}

static void test_compare_uint32(void)
{
	uint32_t a = 0U;
	uint32_t b = 0U;
	uint32_t c = 1U;
	int8_t cmp_eq = vcn_compare_uint32(&a, &b);
	int8_t cmp_lt = vcn_compare_uint32(&b, &c);
	int8_t cmp_gt = vcn_compare_uint32(&c, &a);
	CU_ASSERT(cmp_eq == 0);
	CU_ASSERT(cmp_lt == -1);
	CU_ASSERT(cmp_gt == 1);
}

static void test_compare_uint64(void)
{
	uint64_t a = 0U;
	uint64_t b = 0U;
	uint64_t c = 1U;
	int8_t cmp_eq = vcn_compare_uint64(&a, &b);
	int8_t cmp_lt = vcn_compare_uint64(&b, &c);
	int8_t cmp_gt = vcn_compare_uint64(&c, &a);
	CU_ASSERT(cmp_eq == 0);
	CU_ASSERT(cmp_lt == -1);
	CU_ASSERT(cmp_gt == 1);
}
