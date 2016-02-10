#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <math.h>
#include <string.h>

#include "nb/container_bot/array.h"

#include "complexity_funcs.h"
#include "test.h"
#include "test_driver.h"

#define OUTFILE "complexity_tests.out"
#define LOGFILE "complexity_tests.log"
#define MAX_N_DATA_POINTS (200)
#define MIN(a,b) (((a)<(b))?(a):(b))

typedef struct {
	char* label;
	double (*eval)(double x, double coeff);
} cfunc_t;

enum {
	COMPLEXITY_BIG_O_1 = 0,
	COMPLEXITY_BIG_O_LOGN,
	COMPLEXITY_BIG_O_LOGN2,
	COMPLEXITY_BIG_O_LOGN3,
	COMPLEXITY_BIG_O_N,
	COMPLEXITY_BIG_O_NLOGN,
	COMPLEXITY_BIG_O_NLOGN2,
	COMPLEXITY_BIG_O_NLOGN3,
	COMPLEXITY_BIG_O_N2,
	COMPLEXITY_BIG_O_N2LOGN,
	COMPLEXITY_BIG_O_N2LOGN2,
	COMPLEXITY_BIG_O_N2LOGN3,
	COMPLEXITY_BIG_O_N3,
	COMPLEXITY_BIG_O_N3LOGN,
	COMPLEXITY_BIG_O_N3LOGN2,
	COMPLEXITY_BIG_O_N3LOGN3,
	COMPLEXITY_BIG_O_N4,
	COMPLEXITY_BIG_O_EXPN,
	COMPLEXITY_END
};

static void show_log_opening_test(const char *test_name);
static void show_log(const char *message, int level);
static cfunc_t* estimate_complexity(test_t *test);
static double* get_points(test_t *test, int N_points);
static void show_log_processing_items(int i, int N, int N_items);
static double normalise_points(int N_points, double *points);
static double* add_N_points(test_t *test, int N_points, double *points,
			    int N_add, double normaliser);
static cfunc_t* cfunc_create(void);
static void cfunc_set_data(cfunc_t *cfunc, int complexity_id);
static double fit_coeff(int N_points, double *points, cfunc_t *type);
static double get_error_using_coeff(int N_points, double *points,
				    double coeff, cfunc_t *type);
static void copy_cfunc_t(cfunc_t *dest, cfunc_t *src);
static bool is_accurate(int32_t imin, double error[]);
static void show_complexity_of_test(char *label, const char *bigO_label);
static void cfunc_destroy(cfunc_t *cfunc);
static const char* cfunc_get_label(cfunc_t *cfunc);

int vcn_test_get_driver_type(void)
{
	return VCN_DRIVER_COMPLEXITY_TEST;
}

void vcn_test_do_before_tests(void)
{
	FILE *fp = fopen(OUTFILE, "w");
	fprintf(fp, "ALGORITHMIC COMPLEXITY ESTIMATIONS.\n\n");
	fprintf(fp, "   Big O\tName\n");
	fclose(fp);

	fp = fopen(LOGFILE, "w");
	fprintf(fp, "LOGFILE: ALGORITHMIC COMPLEXITY TESTS.\n\n");
	fclose(fp);
}

bool vcn_test_execute_test(test_t *test)
{
	show_log_opening_test(test->label);
	cfunc_t *bigO = estimate_complexity(test);
	show_complexity_of_test(test->label, cfunc_get_label(bigO));
	cfunc_destroy(bigO);
	show_log("Finish.\n", 0);
	return true;
}

static void show_log_opening_test(const char *test_name)
{
	char log_msg[50];
	sprintf(log_msg, "Testing '%s' ...", test_name);
	show_log(log_msg, 0);
}

static void show_log(const char *message, int level)
{
	FILE *fp = fopen(LOGFILE, "a");
	for (int i = 0; i < level; i++)
		fprintf(fp, "  ");
	fprintf(fp, "%s\n", message);
	fclose(fp);
}

static cfunc_t* estimate_complexity(test_t *test)
{
	int N_points = 10;
	double *points = get_points(test, N_points);
	double normaliser = normalise_points(N_points, points);

	show_log("Estimating complexity...", 1);
	cfunc_t *bigO = cfunc_create();
	double error[COMPLEXITY_END];
	int32_t imin;
	do {
		points = add_N_points(test, N_points, points, 10, normaliser);
		N_points += 10;		
		for (int16_t i = COMPLEXITY_BIG_O_1; i < COMPLEXITY_END; i++) {
			cfunc_set_data(bigO, i);
			double coeff = fit_coeff(N_points, points, bigO);
			error[i] = 1e6;
			if (coeff >= 1.0 && coeff <= 9.0)
				error[i] = get_error_using_coeff(N_points, points, coeff, bigO);
		}
		imin = vcn_array_get_min_id(error, COMPLEXITY_END,
					    sizeof(*error), vcn_compare_double);
	} while (!is_accurate(imin, error) && N_points < MAX_N_DATA_POINTS);

	free(points);
	cfunc_set_data(bigO, imin);
	return bigO;
}

static double* get_points(test_t *test, int N_points)
{
	show_log("Generating data points...", 1);
	double *x = malloc(N_points * 2 * sizeof(*x));
	for (int i = 0; i < N_points; i++) {
		int N_items = (i + 1) * 10;
		show_log_processing_items(i+1, N_points, N_items);
		double seconds = test_get_avg_seconds(test, &N_items);
		x[i * 2] = N_items;
		x[i*2+1] = seconds;
	}
	return x;
}

static void show_log_processing_items(int i, int N, int N_items)
{
	char log_msg[50];
	sprintf(log_msg, "[%i/%i] Running with %i items...", i, N, N_items);
	show_log(log_msg, 2);
}

static double normalise_points(int N_points, double *points)
{
	double mult = 1.0;
	if (N_points > 0) {
		double factor = 5.0; /* Empirically tunned */
		mult = factor * points[0] / points[1];
		for (int i = 0; i < N_points; i++)
			points[i*2+1] *= mult;
	}
	return mult;
}

static double* add_N_points(test_t *test, int N_points, double *points,
			    int N_add, double normaliser)
{
	show_log("Generating more data points to increase accuracy...", 1);
	int N = N_points + N_add;
	double *x = malloc(N * 2 * sizeof(*x));
	memcpy(x, points, 2 * N_points * sizeof(*x));
	free(points);
	for (int i = N_points; i < N; i++) {
		int N_items = (i + 1) * 10;
		show_log_processing_items(i+1, N, N_items);
		double seconds = test_get_avg_seconds(test, &N_items);
		x[i * 2] = N_items;
		x[i*2+1] = seconds * normaliser;
	}
	return x;
}

static inline cfunc_t* cfunc_create(void)
{
	return calloc(1, sizeof(cfunc_t));
}

static void cfunc_set_data(cfunc_t *cfunc, int complexity_id)
{
	switch (complexity_id) {
	case COMPLEXITY_BIG_O_1:
		cfunc->label = "O(1)       ";
		cfunc->eval = eval_1;
		break;
	case COMPLEXITY_BIG_O_LOGN:
		cfunc->label = "O(log(n))  ";
		cfunc->eval = eval_logn;
		break;
	case COMPLEXITY_BIG_O_LOGN2:
		cfunc->label = "O(log(n^2)) ";
		cfunc->eval = eval_logn2;
		break;
	case COMPLEXITY_BIG_O_LOGN3:
		cfunc->label = "O(log(n^3))";
		cfunc->eval = eval_logn3;
		break;
	case COMPLEXITY_BIG_O_N:
		cfunc->label = "O(n)       ";
		cfunc->eval = eval_n;
		break;
	case COMPLEXITY_BIG_O_NLOGN:
		cfunc->label = "O(nlog(n))";
		cfunc->eval = eval_nlogn;
		break;
	case COMPLEXITY_BIG_O_NLOGN2:
		cfunc->label = "O(nlog(n^2))";
		cfunc->eval = eval_nlogn2;
		break;
	case COMPLEXITY_BIG_O_NLOGN3:
		cfunc->label = "O(nlog(n^3))";
		cfunc->eval = eval_nlogn3;
		break;
	case COMPLEXITY_BIG_O_N2:
		cfunc->label = "O(n^2)      ";
		cfunc->eval = eval_n2;
		break;
	case COMPLEXITY_BIG_O_N2LOGN:
		cfunc->label = "O(n^2 log(n))";
		cfunc->eval = eval_n2logn;
		break;
	case COMPLEXITY_BIG_O_N2LOGN2:
		cfunc->label = "O(n^2 log(n^2))";
		cfunc->eval = eval_n2logn2;
		break;
	case COMPLEXITY_BIG_O_N2LOGN3:
		cfunc->label = "O(n^2 log(n^3))";
		cfunc->eval = eval_n2logn3;
		break;
	case COMPLEXITY_BIG_O_N3:
		cfunc->label = "O(n^3)      ";
		cfunc->eval = eval_n3;
		break;
	case COMPLEXITY_BIG_O_N3LOGN:
		cfunc->label = "O(n^3 log(n))";
		cfunc->eval = eval_n3logn;
		break;
	case COMPLEXITY_BIG_O_N3LOGN2:
		cfunc->label = "O(n^3 log(n^2))";
		cfunc->eval = eval_n3logn2;
		break;
	case COMPLEXITY_BIG_O_N3LOGN3:
		cfunc->label = "O(n^3 log(n^3))";
		cfunc->eval = eval_n3logn3;
		break;
	case COMPLEXITY_BIG_O_N4:
		cfunc->label = "O(n^4)      ";
		cfunc->eval = eval_n4;
		break;
	case COMPLEXITY_BIG_O_EXPN:
		cfunc->label = "O(exp(n))   ";
		cfunc->eval = eval_expn;
		break;
	}
}

static double fit_coeff(int N_points, double *points, cfunc_t *type)
{
	double sum_fi2 = 0;
	double sum_fiyi = 0;
	for (int i = 0; i < N_points; i++) {
		double fi = type->eval(points[i*2], 1);
		sum_fi2 += fi * fi;
		sum_fiyi += fi * points[i*2+1];
	}
	return sum_fiyi / sum_fi2;
}

static double get_error_using_coeff(int N_points, double *points,
				    double coeff, cfunc_t *type)
{
	double error_sum = 0;
	for (int i = 0; i < N_points; i++) {
		double eval = type->eval(points[i * 2], coeff);
		double error = eval - points[i*2+1];
		error_sum += error * error;
	}
	return sqrt(error_sum);
}

static inline void copy_cfunc_t(cfunc_t *dest, cfunc_t *src)
{
	memcpy(dest, src, sizeof(*src));
}

static bool is_accurate(int32_t imin, double error[])
{
	bool is_accur = true;
	for (int16_t i = COMPLEXITY_BIG_O_1; i < COMPLEXITY_END; i++) {
		if (i != imin) {
			if (fabs(error[imin] - error[i]) < 1.0) {
				is_accur = false;
				break;
			}
		}
	}
	return is_accur;
}

static inline void show_complexity_of_test(char *label, const char *bigO_label)
{
	FILE *fp = fopen(OUTFILE, "a");
	fprintf(fp, " > %s\t%s\n", bigO_label, label);
	fclose(fp);
}

static inline const char* cfunc_get_label(cfunc_t *cfunc)
{
	return cfunc->label;
}

static inline void cfunc_destroy(cfunc_t *cfunc)
{
	free(cfunc);
}
