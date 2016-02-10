#include <stdlib.h>

#include "nb/container_bot/container.h"

#include "test.h"
#include "test_add.h"

static void test_add(vcn_container_t *tests,
		     void* (*init)(),
		     bool (*run)(),
		     void (*clear)(void*),
		     char *label);

static void* null_init();
static void null_clear(void* data);

inline void vcn_test_add(void *tests_ptr,
			 bool (*run)(void),
			 char *label)
{
	test_add(tests_ptr, null_init, run, null_clear, label);
}

static void test_add(vcn_container_t *tests,
		     void* (*init)(),
		     bool (*run)(),
		     void (*clear)(void*),
		     char *label)
{
	test_t *test = test_create();
	test->label = label;
	test->init = init;
	test->run = run;
	test->clear = clear;
	vcn_container_insert(tests, test);
}

static inline void* null_init()
{
	return NULL;
}

static inline void null_clear(void* data)
{
	; /* Null statement */
}

inline void vcn_test_add_decomposed(void *tests_ptr,
				    void* (*init)(),
				    bool (*run)(void*),
				    void (*clear)(void*),
				    char *label)
{
	test_add(tests_ptr, init, run, clear, label);
}
