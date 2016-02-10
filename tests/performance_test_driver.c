#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "test.h"
#include "test_driver.h"

#define OUTFILE "performance_tests.out"

static void show_performance_of_test(char *label, double sec);

inline int vcn_test_get_driver_type(void)
{
	return NB_DRIVER_PERFORMANCE_TEST;
}

void vcn_test_do_before_tests(void)
{
	FILE *fp = fopen(OUTFILE, "w");
	fprintf(fp, "   Seconds\tName\n");
	fclose(fp);
}

bool vcn_test_execute_test(test_t *test)
{
  double time = test_get_avg_seconds(test, NULL);
	show_performance_of_test(test->label, time);
	return true;
}

static inline void show_performance_of_test(char *label, double sec)
{
	FILE *fp = fopen(OUTFILE, "a");
	fprintf(fp, " > %.3e\t%s\n", sec, label);
	fclose(fp);
}
