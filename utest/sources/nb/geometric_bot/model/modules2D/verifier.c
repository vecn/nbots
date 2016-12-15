#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "cunit/Basic.h"

#include  "nb/math_bot.h"
#include  "nb/geometric_bot/mesh/mesh2D.h"
#include  "nb/geometric_bot/model/modules2D/verifier.h"

#define INPUTS_DIR "../utest/sources/nb/geometric_bot/model/modules2D/verifier_inputs"

static int suite_init(void);
static int suite_clean(void);

static void test_verify_consistence(void);

void cunit_nb_geometric_bot_model2D_verifier(void)
{
	CU_pSuite suite = CU_add_suite("nb/geometric_bot/model/modules2D/verifier.c",
				       suite_init, suite_clean);
	CU_add_test(suite, "verify_consistence() of square-donut",
		    test_verify_consistence);
}

static int suite_init(void)
{
	return 0;
}

static int suite_clean(void)
{
	return 0;
}

static void test_verify_consistence(void)
{
	char input_name[256];
	sprintf(input_name, "%s/square_donut.psl", INPUTS_DIR);
	nb_model_t *model = nb_model_load(input_name);
	bool is_ok = (0 == nb_model_verify_consistence(model, NULL));
	nb_model_destroy(model);
	CU_ASSERT(is_ok);
}
