#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "test_library.h"
#include "test_add.h"

#include  "nb/math_bot.h"
#include  "nb/geometric_bot/mesh/mesh2D.h"
#include  "nb/geometric_bot/model/modules2D/verifier.h"

#define INPUTS_DIR "../tests/nb/geometric_bot/model/modules2D/verifier_UT_inputs"

static bool check_verify_consistence(void);

inline int vcn_test_get_driver_id(void)
{
	return NB_DRIVER_UNIT_TEST;
}

void vcn_test_load_tests(void *tests_ptr)
{
	vcn_test_add(tests_ptr, check_verify_consistence,
		     "Check verify_consistence() of square-donut");
}

static bool check_verify_consistence(void)
{
	char input_name[256];
	sprintf(input_name, "%s/square_donut.psl", INPUTS_DIR);
	vcn_model_t *model = vcn_model_load(input_name);
	bool is_ok = (0 == vcn_model_verify_consistence(model, NULL));
	vcn_model_destroy(model);
	return is_ok;
}
