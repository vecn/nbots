#include "cunit/CUnit.h"
#include "cunit/Basic.h"
#include "nbots.h"

#include "nb/container_bot_utests.h"
#include "nb/geometric_bot_utests.h"
#include "nb/image_bot_utests.h"
#include "nb/pde_bot_utests.h"
#include "nb/graphics_bot_utests.h"


int main(int argc, char *argv[])
{
	CU_initialize_registry();

	cunit_suites_nb_container_bot();
	cunit_suites_nb_geometric_bot();
	cunit_suites_nb_image_bot();
	cunit_suites_nb_graphics_bot();
	cunit_suites_nb_pde_bot();
	
	CU_basic_set_mode(CU_BRM_VERBOSE);
	CU_basic_run_tests();
	CU_cleanup_registry();
	return CU_get_error();
}
