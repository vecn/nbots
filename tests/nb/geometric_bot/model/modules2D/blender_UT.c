#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "test_library.h"
#include "test_add.h"

#include  "nb/math_bot.h"
#include  "nb/geometric_bot/mesh/mesh2D.h"

#define INPUTS_DIR "../tests/nb/geometric_bot/model/modules2D/blender_UT_inputs"

#include  "nb/geometric_bot/model/modules2D/exporter_cairo.h" /* TEMPORAL */
static int TEMPORAL_ = 0; /* TEMPORAL */		      /* TEMPORAL */
static void TEMPORAL(const vcn_model_t *const model)	      /* TEMPORAL */
{							      /* TEMPORAL */
	char label[100];				      /* TEMPORAL */
	sprintf(label, "TEMP_M_%02i.png", TEMPORAL_++);	      /* TEMPORAL */
	vcn_model_export_png(model, label, 1000, 800, false);
	//vcn_mesh_t *mesh = vcn_mesh_create();
	//vcn_mesh_generate_from_model(mesh, model);
	//vcn_mesh_save_png(mesh, label, 1000, 800);  /* TEMPORAL */
	//vcn_mesh_destroy(mesh);
}                                                             /* TEMPORAL */

static bool check_get_combination_squares(void);
static bool check_get_intersection_squares(void);
static bool check_get_union_squares(void);
static bool check_get_substractionA_squares(void);
static bool check_get_substractionB_squares(void);
static bool check_get_combination_cross(void);
static bool check_get_intersection_cross(void);
static bool check_get_union_cross(void);
static bool check_get_substractionA_cross(void);
static bool check_get_substractionB_cross(void);

inline int vcn_test_get_driver_id(void)
{
	return NB_DRIVER_UNIT_TEST;
}

void vcn_test_load_tests(void *tests_ptr)
{
	//vcn_test_add(tests_ptr, check_get_combination_squares,
	//	     "Check get_combination() of squares");
	//vcn_test_add(tests_ptr, check_get_intersection_squares,
	//	     "Check get_intersection() of squares");
	//vcn_test_add(tests_ptr, check_get_union_squares,
	//	     "Check get_union() of squares");
	//vcn_test_add(tests_ptr, check_get_substractionA_squares,
	//	     "Check get_substraction() of squares A");
	//vcn_test_add(tests_ptr, check_get_substractionB_squares,
	//	     "Check get_substraction() of squares B");
	vcn_test_add(tests_ptr, check_get_combination_cross,
		     "Check get_combination() of cross");
	vcn_test_add(tests_ptr, check_get_intersection_cross,
		     "Check get_intersection() of cross");
	vcn_test_add(tests_ptr, check_get_union_cross,
		     "Check get_union() of cross");
	vcn_test_add(tests_ptr, check_get_substractionA_cross,
		     "Check get_substraction() of cross A");
	vcn_test_add(tests_ptr, check_get_substractionB_cross,
		     "Check get_substraction() of cross B");
}

static bool check_get_combination_squares(void)
{
	char input_name[256];
	sprintf(input_name, "%s/combining_squares_1.psl", INPUTS_DIR);
	vcn_model_t *model1 = vcn_model_load(input_name);
	sprintf(input_name, "%s/combining_squares_2.psl", INPUTS_DIR);
	vcn_model_t *model2 = vcn_model_load(input_name);
	vcn_model_t *model = vcn_model_get_combination(model1, model2, 0);
	TEMPORAL(model); /* TEMPORAL */
	vcn_model_destroy(model1);
	vcn_model_destroy(model2);
	vcn_model_destroy(model);
	return false;
}

static bool check_get_intersection_squares(void)
{
	char input_name[256];
	sprintf(input_name, "%s/combining_squares_1.psl", INPUTS_DIR);
	vcn_model_t *model1 = vcn_model_load(input_name);
	sprintf(input_name, "%s/combining_squares_2.psl", INPUTS_DIR);
	vcn_model_t *model2 = vcn_model_load(input_name);
	vcn_model_t *model = vcn_model_get_intersection(model1, model2, 0);
	TEMPORAL(model); /* TEMPORAL */
	vcn_model_destroy(model1);
	vcn_model_destroy(model2);
	vcn_model_destroy(model);
	return false;
}

static bool check_get_union_squares(void)
{
	char input_name[256];
	sprintf(input_name, "%s/combining_squares_1.psl", INPUTS_DIR);
	vcn_model_t *model1 = vcn_model_load(input_name);
	sprintf(input_name, "%s/combining_squares_2.psl", INPUTS_DIR);
	vcn_model_t *model2 = vcn_model_load(input_name);
	vcn_model_t *model = vcn_model_get_union(model1, model2, 0);
	TEMPORAL(model); /* TEMPORAL */
	vcn_model_destroy(model1);
	vcn_model_destroy(model2);
	vcn_model_destroy(model);
	return false;
}

static bool check_get_substractionA_squares(void)
{
	char input_name[256];
	sprintf(input_name, "%s/combining_squares_1.psl", INPUTS_DIR);
	vcn_model_t *model1 = vcn_model_load(input_name);
	sprintf(input_name, "%s/combining_squares_2.psl", INPUTS_DIR);
	vcn_model_t *model2 = vcn_model_load(input_name);
	vcn_model_t *model = vcn_model_get_substraction(model1, model2, 0);
	TEMPORAL(model); /* TEMPORAL */
	vcn_model_destroy(model1);
	vcn_model_destroy(model2);
	vcn_model_destroy(model);
	return false;
}

static bool check_get_substractionB_squares(void)
{
	char input_name[256];
	sprintf(input_name, "%s/combining_squares_2.psl", INPUTS_DIR);
	vcn_model_t *model1 = vcn_model_load(input_name);
	sprintf(input_name, "%s/combining_squares_1.psl", INPUTS_DIR);
	vcn_model_t *model2 = vcn_model_load(input_name);
	vcn_model_t *model = vcn_model_get_substraction(model1, model2, 0);
	TEMPORAL(model); /* TEMPORAL */
	vcn_model_destroy(model1);
	vcn_model_destroy(model2);
	vcn_model_destroy(model);
	return false;
}

static bool check_get_combination_cross(void)
{
	char input_name[256];
	sprintf(input_name, "%s/combining_cross_1.psl", INPUTS_DIR);
	vcn_model_t *model1 = vcn_model_load(input_name);
	sprintf(input_name, "%s/combining_cross_2.psl", INPUTS_DIR);
	vcn_model_t *model2 = vcn_model_load(input_name);
	vcn_model_t *model = vcn_model_get_combination(model1, model2, 0);
	TEMPORAL(model); /* TEMPORAL */
	vcn_model_destroy(model1);
	vcn_model_destroy(model2);
	vcn_model_destroy(model);
	return false;
}

static bool check_get_intersection_cross(void)
{
	char input_name[256];
	sprintf(input_name, "%s/combining_cross_1.psl", INPUTS_DIR);
	vcn_model_t *model1 = vcn_model_load(input_name);
	sprintf(input_name, "%s/combining_cross_2.psl", INPUTS_DIR);
	vcn_model_t *model2 = vcn_model_load(input_name);
	vcn_model_t *model = vcn_model_get_intersection(model1, model2, 0);
	TEMPORAL(model); /* TEMPORAL */
	vcn_model_destroy(model1);
	vcn_model_destroy(model2);
	vcn_model_destroy(model);
	return false;
}

static bool check_get_union_cross(void)
{
	char input_name[256];
	sprintf(input_name, "%s/combining_cross_1.psl", INPUTS_DIR);
	vcn_model_t *model1 = vcn_model_load(input_name);
	sprintf(input_name, "%s/combining_cross_2.psl", INPUTS_DIR);
	vcn_model_t *model2 = vcn_model_load(input_name);
	vcn_model_t *model = vcn_model_get_union(model1, model2, 0);
	TEMPORAL(model); /* TEMPORAL */
	vcn_model_destroy(model1);
	vcn_model_destroy(model2);
	vcn_model_destroy(model);
	return false;
}

static bool check_get_substractionA_cross(void)
{
	char input_name[256];
	sprintf(input_name, "%s/combining_cross_1.psl", INPUTS_DIR);
	vcn_model_t *model1 = vcn_model_load(input_name);
	sprintf(input_name, "%s/combining_cross_2.psl", INPUTS_DIR);
	vcn_model_t *model2 = vcn_model_load(input_name);
	vcn_model_t *model = vcn_model_get_substraction(model1, model2, 0);
	TEMPORAL(model); /* TEMPORAL */
	vcn_model_destroy(model1);
	vcn_model_destroy(model2);
	vcn_model_destroy(model);
	return false;
}

static bool check_get_substractionB_cross(void)
{
	char input_name[256];
	sprintf(input_name, "%s/combining_cross_2.psl", INPUTS_DIR);
	vcn_model_t *model1 = vcn_model_load(input_name);
	sprintf(input_name, "%s/combining_cross_1.psl", INPUTS_DIR);
	vcn_model_t *model2 = vcn_model_load(input_name);
	vcn_model_t *model = vcn_model_get_substraction(model1, model2, 0);
	TEMPORAL(model); /* TEMPORAL */
	vcn_model_destroy(model1);
	vcn_model_destroy(model2);
	vcn_model_destroy(model);
	return false;
}
