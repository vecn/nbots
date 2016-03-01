#define _XOPEN_SOURCE /* To compile sigaction() */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>

#include <dlfcn.h>  /* POSIX */
#include <signal.h> /* POSIX */

#include "nbots.h"

#include "test.h"
#include "test_driver.h"

enum {
	EXCEPTION_NULL_LIBRARY,
	EXCEPTION_NULL_DRIVER,
	EXCEPTION_WRONG_DRIVER,
	EXCEPTION_NULL_LOADER
};

void *signal_data;

static void define_signals_behaviour();
static void set_signal_handler(int signal, void (*handler)(int));
static void handle_signal(int signal);
static char* get_signal_str(int signal);
static bool try_to_process_library(const char *library_name);
static bool process_library_or_throw_exception(vcn_exception_t *exception,
					       const char *library_name);
static void* load_library_or_throw_exception(vcn_exception_t *exception,
					     const char *library_name);
static void throw_exception_if_lib_driver_is_wrong(vcn_exception_t *exception,
						   void *lib);
static void* get_loader_or_throw_exception(vcn_exception_t *exception,
					   void *lib);
static void destroy_tests(void* tests);
static bool process_tests(nb_container_t *tests);
static int execute_tests(int N_tests, void **tests);
static void show_running_msg(int i, int N, const char *label);
static void clear_stream(FILE *stream);
static int try_to_execute_test_and_catch_signals(test_t *test);
static int execute_test_and_count_if_success(test_t *test);
static void show_msg_after_tests(int *data);
static void handle_exception(vcn_exception_t *exception);
static void close_library(void *lib);

int main(int argc, char *argv[])
{
	int exit_status;
	if (argc < 2) {
		const char *error_label =
			"Provide the tests-library name as the 1st argument";
		printf("ERROR: %s.\n", error_label);
		exit_status = 1;
	} else {
		define_signals_behaviour();
		bool all_pass = try_to_process_library(argv[1]);
		exit_status = (all_pass)?0:2;
	}
	return exit_status;
}

static void define_signals_behaviour()
{
	set_signal_handler(SIGSEGV, handle_signal);
	set_signal_handler(SIGTERM, handle_signal);
	set_signal_handler(SIGINT, handle_signal);
	set_signal_handler(SIGILL, handle_signal);
	set_signal_handler(SIGABRT, handle_signal);
	set_signal_handler(SIGFPE, handle_signal);
}

static void set_signal_handler(int signal, void (*handler)(int))
{
	struct sigaction act;
	memset(&act, 0, sizeof(act));
	act.sa_handler = handler;
	if (sigaction(signal, &act, NULL) == -1)
		printf("Can not handle signal %i.\n", signal);
}

static void handle_signal(int signal)
{
	vcn_exception_t *exception = signal_data;
	vcn_exception_throw(exception, signal, get_signal_str(signal));
}

static char* get_signal_str(int signal)
{
	char *str;
	switch (signal) {
	case SIGSEGV:
		str = "Segmentation fault";
		break;
	case SIGTERM:
		str = "Termination request";
		break;
	case SIGINT:
		str = "External interrupt";
		break;
	case SIGILL:
		str = "Invalid instruction";
		break;
	case SIGABRT:
		str = "Abnormal termination, e.g. abort()";
		break;
	case SIGFPE:
		str = "Floating point exception";
		break;
	default:
		str = "Unknown signal";
		break;
	}
	return str;
}

static bool try_to_process_library(const char *library_name)
{
	bool all_pass = false;
	vcn_exception_t *exception = vcn_exception_create();
	if(vcn_exception_try(exception))
		all_pass = process_library_or_throw_exception(exception,
							      library_name);
	else
		handle_exception(exception);
	vcn_exception_clear_alloc(exception);
	vcn_exception_destroy(exception);
	return all_pass;
}

static bool process_library_or_throw_exception(vcn_exception_t *exception,
					       const char *library_name)
{
	void *lib = load_library_or_throw_exception(exception, library_name);
	vcn_exception_set_alloc(exception, lib, close_library);

	throw_exception_if_lib_driver_is_wrong(exception, lib);

	void (*load_tests)(nb_container_t*);
	load_tests = get_loader_or_throw_exception(exception, lib);

	nb_container_t *tests = nb_container_create(NB_QUEUE);
	vcn_exception_set_alloc(exception, tests, destroy_tests);
	load_tests(tests);

	vcn_test_do_before_tests();

	bool all_pass = process_tests(tests);
	return all_pass;
}

static void* load_library_or_throw_exception(vcn_exception_t *exception,
					     const char *library_name)
/* POSIX systems */
{
	void *lib = dlopen(library_name, RTLD_NOW);
	if (NULL == lib)
		vcn_exception_throw(exception, EXCEPTION_NULL_LIBRARY,
				    dlerror());
	return lib;
}

static void close_library(void *lib)
/* POSIX systems */
{
	dlclose(lib);
}

static void throw_exception_if_lib_driver_is_wrong(vcn_exception_t *exception,
						   void *lib)
/* POSIX systems */
{  
	void* get_driver_ptr = dlsym(lib, "vcn_test_get_driver_id");
	if (NULL == get_driver_ptr)
		vcn_exception_throw(exception, EXCEPTION_NULL_DRIVER, dlerror());

	int (*get_driver_to_process_library)() = get_driver_ptr;
	if (vcn_test_get_driver_type() != get_driver_to_process_library())
		vcn_exception_throw(exception, EXCEPTION_WRONG_DRIVER,
				    "wrong driver");
}

static void* get_loader_or_throw_exception(vcn_exception_t *exception, void *lib)
/* POSIX systems */
{
	void* loader_ptr = dlsym(lib, "vcn_test_load_tests");
	if (NULL == loader_ptr)
		vcn_exception_throw(exception, EXCEPTION_NULL_LOADER, dlerror());
	return loader_ptr;
}

static void destroy_tests(void* tests)
{
	nb_container_destroy(tests);
}

static bool process_tests(nb_container_t *tests)
{
	int N_tests = nb_container_get_length(tests);
	int N_success = 0;
	if (N_tests > 0) {
		void** tests_array = malloc(N_tests * sizeof(*tests_array));
		nb_container_copy_to_array(tests, tests_array);
		N_success += execute_tests(N_tests, tests_array);
		free(tests_array);
	}
	int data[2];
	data[0] = N_success;
	data[1] = N_tests;
	show_msg_after_tests(data);
	return (N_success == N_tests);
}

static int execute_tests(int N_tests, void **tests)
{
	int N_success = 0;
	uint32_t *tests_permuted = malloc(N_tests * sizeof(*tests_permuted));
	for (int i = 0; i < N_tests; i++)
		tests_permuted[i] = i;
	vcn_statistics_random_permutation(N_tests, tests_permuted,
					  sizeof(*tests_permuted));
	for (int i = 0; i < N_tests; i++) {
		int id = tests_permuted[i];
		test_t *test = tests[id];
		show_running_msg(i + 1, N_tests, test->label);
		N_success += try_to_execute_test_and_catch_signals(test);
		clear_stream(stdout);
	}
	free(tests_permuted);
	return N_success;
}

static void show_running_msg(int i, int N, const char *label)
{
	printf("\r > [%i/%i] Running Test: '%s' ...", i, N, label);
	fflush(stdout);

}

static void clear_stream(FILE *stream)
{
	char clear_str[100];
	memset(clear_str, ' ', 100);
	fprintf(stream, "\r%s", clear_str);
	fflush(stream);
}

static int try_to_execute_test_and_catch_signals(test_t *test)
{
	vcn_exception_t* exception_by_signal = vcn_exception_create();
	signal_data = exception_by_signal;
	int count = 0;
	if (vcn_exception_try(exception_by_signal)) {
		count = execute_test_and_count_if_success(test);
	} else {
		printf("\r > Test fails: '%s'  (%s).        \n", 
		       test->label,
		       vcn_exception_get_info(exception_by_signal));
		exit(1);
	}
	vcn_exception_destroy(exception_by_signal);
	return count;
}

static int execute_test_and_count_if_success(test_t *test)
{
	bool test_pass = vcn_test_execute_test(test);
	bool test_fails = !test_pass;
	if (test_fails)
		printf("\r > Test fails: '%s'.              \n", test->label);
	return (test_pass)?1:0;
}

static void show_msg_after_tests(int *data)
{
	int N_success = data[0];
	int N_total = data[1];
	if (0 == N_total)
		printf("\r[%i/%i]--------> Zero tests\n", N_success, N_total);
	else if (N_success == N_total)
		printf("\r[%i/%i] Success\n", N_success, N_total);
	else
		printf("\r[%i/%i]--------> Something wrong\n", 
		       N_success, N_total);
}

static void handle_exception(vcn_exception_t *exception)
{
	switch (vcn_exception_get_id(exception)) {
	case EXCEPTION_NULL_LIBRARY:
		printf("Failed to open tests-library (%s).\n",
		       vcn_exception_get_info(exception));
		break;
	case EXCEPTION_NULL_DRIVER:
		
		printf("Failed to locate vcn_test_get_driver_id() (%s).\n",
		       vcn_exception_get_info(exception));
		break;
	case EXCEPTION_WRONG_DRIVER:
		printf("The library requires another test driver (%s).\n",
		       vcn_exception_get_info(exception));
		break;
	case EXCEPTION_NULL_LOADER:
		printf("Failed to locate vcn_test_load_tests(). (%s)\n",
		       vcn_exception_get_info(exception));
		break;
	default:
		printf("Unhandled exception occurs\n");
		break;
	}
}
