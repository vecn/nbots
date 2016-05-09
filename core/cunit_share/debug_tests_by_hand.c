/* Compile from build dir:
  export LD_LIBRARY_PATH=libs/nbots/shared/debug:libs/nbots_cairo/shared/debug
  gcc -std=c99 -g -I ../include/ -L libs/nbots/shared/debug -L libs/nbots_cairo/shared/debug ../cunit_share/debug_tests_by_hand.c -o debug_tests_by_hand -lm -lnbots -lnbots_cairo -lcunit
 */

#include "nbots.h"

#include "../utest/nb/container_bot_utests.h"
#include "../utest/nb/geometric_bot_utests.h"
#include "../utest/nb/pde_bot_utests.h"

#include <CUnit/Basic.h>

int main()
{
	cunit_suites_nb_container_bot();
	cunit_suites_nb_geometric_bot();
	cunit_suites_nb_pde_bot();

	if (CUE_SUCCESS == CU_initialize_registry()) {
		/* Run all tests using the CUnit Basic interface */
		CU_basic_set_mode(CU_BRM_VERBOSE);
		CU_basic_run_tests();
		CU_cleanup_registry();
	} else {
		printf("ERROR: Initializing CUnit registry...");
	}
	return 0;
}
