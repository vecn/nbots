#ifndef __VCN_TEST_DRIVER_H__
#define __VCN_TEST_DRIVER_H__

#include <stdbool.h>
#include "test.h"
#include "test_driver_types.h"

int vcn_test_get_driver_type(void);
void vcn_test_do_before_tests(void);
bool vcn_test_execute_test(test_t *test);

#endif
