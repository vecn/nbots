#define CONTAINER_ID VCN_CONTAINER_SORTED
#define N_ITEMS 5000
#include "iterator_TEMPLATE_UT.c"

inline int vcn_test_get_driver_id(void)
{
	return TEMPLATE_get_driver_id();
}

void vcn_test_load_tests(void *tests_ptr)
{
	TEMPLATE_load_tests(tests_ptr);
}
