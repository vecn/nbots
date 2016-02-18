#include <stdlib.h>
#include "test.h"

#include <sys/time.h>     /* Linux kernel */

static double get_total_seconds(const test_t *test, int N_runs, void *init_data);
static int get_N_runs(double time);
static double get_seconds(struct timeval start, struct timeval end);

test_t* test_create(void)
{
	return calloc(1, sizeof(test_t));
}

double test_get_avg_seconds(const test_t *test, void *init_data)
{
	int N_runs = 5;
	double time = get_total_seconds(test, N_runs, init_data);
	double avg_time = time / N_runs;
	if (avg_time < 1e-2) {
	  N_runs = get_N_runs(avg_time);
	  time = get_total_seconds(test, N_runs, init_data);
	}
	return time / N_runs;
}

static double get_total_seconds(const test_t *test, int N_runs, void *init_data)
{
	struct timeval start, end;
	double running_time = 0.0;
	for (int i = 0; i < N_runs; i++) {
		void *test_data = test->init(init_data);

		gettimeofday(&start, NULL);
		test->run(test_data);
		gettimeofday(&end, NULL);

		test->clear(test_data);

		running_time += get_seconds(start, end);
	}
	return running_time;
}

static int get_N_runs(double time)
{
	int N_runs = 50;
	if (time < 1e-3)
		N_runs = 500;
	else if (time < 1e-4)
		N_runs = 10000;
	return N_runs;
}

static double get_seconds(struct timeval start, struct timeval end){
	return (double)((end.tv_sec*1000000 + end.tv_usec) -
			(start.tv_sec*1000000 + start.tv_usec)) / 1e6;
}
