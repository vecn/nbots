#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "test.h"
#include "test_driver.h"

int vcn_test_get_driver_type(void)
{
	return VCN_DRIVER_UNIT_TEST;
}

void vcn_test_do_before_tests(void)
{
	; /* Null statement */
}

bool vcn_test_execute_test(test_t *test)
{
	return test->run(NULL);
}
