#ifndef __NB_TEST_LIBRARY_ADD_H__
#define __NB_TEST_LIBRARY_ADD_H__

#include <stdbool.h>

void vcn_test_add(void *tests_ptr,
		  bool (*run)(void),
		  char *label);

void vcn_test_add_decomposed(void *tests_ptr,
			     void* (*init)(),
			     bool (*run)(void*),
			     void (*clear)(void*),
			     char *label);

#endif
