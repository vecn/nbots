#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "test_library.h"
#include "test_add.h"

#include  "nb/math_bot.h"
#include  "nb/geometric_bot/mesh/mesh2D.h"
#include  "nb/geometric_bot/model/modules2D/clipper.h"

#define INPUTS_DIR "tests/nb/geometric_bot/model/modules2D/clipper_UT_inputs"
#define OUTPUT_DIR "build"

#include  "nb/geometric_bot/model/modules2D/exporter_cairo.h" /* TEMPORAL */
#include  "nb/geometric_bot/mesh/modules2D/exporter_cairo.h" /* TEMPORAL */
static int TEMPORAL_ = 0; /* TEMPORAL */		      /* TEMPORAL */
static void TEMPORAL(const vcn_model_t *const model)	      /* TEMPORAL */
{							      /* TEMPORAL */
	char label[100];				      /* TEMPORAL */
	sprintf(label, "%s/TEMP_MODEL_%02i.png",
		OUTPUT_DIR, TEMPORAL_++);                    /* TEMPORAL */
	vcn_model_export_png(model, label, 1000, 800, false);
	vcn_mesh_t *mesh = vcn_mesh_create();
	vcn_mesh_generate_from_model(mesh, model);
	sprintf(label, "%s/TEMP_MODEL_%02i_MSH.png", OUTPUT_DIR,
		TEMPORAL_++);	                             /* TEMPORAL */
	vcn_mesh_save_png(mesh, label, 1000, 800);  /* TEMPORAL */
	vcn_mesh_destroy(mesh);
}                                                             /* TEMPORAL */

static bool check_get_combination_squares(void);
static bool check_get_intersection_squares(void);
static bool check_get_union_squares(void);
static bool check_get_substractionA_squares(void);
static bool check_get_substractionB_squares(void);
static bool check_get_difference_squares(void);
static bool check_get_combination_cross(void);
static bool check_get_intersection_cross(void);
static bool check_get_union_cross(void);
static bool check_get_substractionA_cross(void);
static bool check_get_substractionB_cross(void);
static bool check_get_difference_cross(void);
static bool check_get_combination_c(void);
static bool check_get_intersection_c(void);
static bool check_get_union_c(void);
static bool check_get_substractionA_c(void);
static bool check_get_substractionB_c(void);
static bool check_get_difference_c(void);
static bool check_get_combination_b(void);
static bool check_get_intersection_b(void);
static bool check_get_union_b(void);
static bool check_get_substractionA_b(void);
static bool check_get_substractionB_b(void);
static bool check_get_difference_b(void);
static bool check_get_combination_g(void);
static bool check_get_intersection_g(void);
static bool check_get_union_g(void);
static bool check_get_substractionA_g(void);
static bool check_get_substractionB_g(void);
static bool check_get_difference_g(void);
static bool check_get_combination_h(void);
static bool check_get_intersection_h(void);
static bool check_get_union_h(void);
static bool check_get_substractionA_h(void);
static bool check_get_substractionB_h(void);
static bool check_get_difference_h(void);

inline int vcn_test_get_driver_id(void)
{
	return NB_DRIVER_UNIT_TEST;
}

void vcn_test_load_tests(void *tests_ptr)
{
	vcn_test_add(tests_ptr, check_get_combination_squares,
		     "Check get_combination() of squares");
	vcn_test_add(tests_ptr, check_get_intersection_squares,
		     "Check get_intersection() of squares");
	vcn_test_add(tests_ptr, check_get_union_squares,
		     "Check get_union() of squares");
	vcn_test_add(tests_ptr, check_get_substractionA_squares,
		     "Check get_substraction() of squares A");
	vcn_test_add(tests_ptr, check_get_substractionB_squares,
		     "Check get_substraction() of squares B");
	vcn_test_add(tests_ptr, check_get_difference_squares,
		     "Check get_difference() of squares B");
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
	vcn_test_add(tests_ptr, check_get_difference_cross,
		     "Check get_difference() of cross");
	vcn_test_add(tests_ptr, check_get_combination_c,
		     "Check get_combination() of c");
	vcn_test_add(tests_ptr, check_get_intersection_c,
		     "Check get_intersection() of c");
	vcn_test_add(tests_ptr, check_get_union_c,
		     "Check get_union() of c");
	vcn_test_add(tests_ptr, check_get_substractionA_c,
		     "Check get_substraction() of c A");
	vcn_test_add(tests_ptr, check_get_substractionB_c,
		     "Check get_substraction() of c B");
	vcn_test_add(tests_ptr, check_get_difference_c,
		     "Check get_difference() of c");
	vcn_test_add(tests_ptr, check_get_combination_b,
		     "Check get_combination() of b");
	vcn_test_add(tests_ptr, check_get_intersection_b,
		     "Check get_intersection() of b");
	vcn_test_add(tests_ptr, check_get_union_b,
		     "Check get_union() of b");
	vcn_test_add(tests_ptr, check_get_substractionA_b,
		     "Check get_substraction() of b A");
	vcn_test_add(tests_ptr, check_get_substractionB_b,
		     "Check get_substraction() of b B");
	vcn_test_add(tests_ptr, check_get_difference_b,
		     "Check get_difference() of b");
	return;/* TEMPORAL */

	vcn_test_add(tests_ptr, check_get_combination_g,
		     "Check get_combination() of g");
	vcn_test_add(tests_ptr, check_get_intersection_g,
		     "Check get_intersection() of g");
	vcn_test_add(tests_ptr, check_get_union_g,
		     "Check get_union() of g");
	vcn_test_add(tests_ptr, check_get_substractionA_g,
		     "Check get_substraction() of g A");
	vcn_test_add(tests_ptr, check_get_substractionB_g,
		     "Check get_substraction() of g B");
	vcn_test_add(tests_ptr, check_get_difference_g,
		     "Check get_difference() of g");
	vcn_test_add(tests_ptr, check_get_combination_h,
		     "Check get_combination() of h");
	vcn_test_add(tests_ptr, check_get_intersection_h,
		     "Check get_intersection() of h");
	vcn_test_add(tests_ptr, check_get_union_h,
		     "Check get_union() of h");
	vcn_test_add(tests_ptr, check_get_substractionA_h,
		     "Check get_substraction() of h A");
	vcn_test_add(tests_ptr, check_get_substractionB_h,
		     "Check get_substraction() of h B");
	vcn_test_add(tests_ptr, check_get_difference_h,
		     "Check get_difference() of h");
}

static bool check_get_combination_squares(void)
{
	char input_name[256];
	sprintf(input_name, "%s/combining_squares_1.psl", INPUTS_DIR);
	vcn_model_t *model1 = vcn_model_load(input_name);
	sprintf(input_name, "%s/combining_squares_2.psl", INPUTS_DIR);
	vcn_model_t *model2 = vcn_model_load(input_name);
	vcn_model_t *model = vcn_model_create();
	vcn_model_get_combination(model, model1, model2, 0);
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
	vcn_model_t *model = vcn_model_create();
	vcn_model_get_intersection(model, model1, model2, 0);
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
	vcn_model_t *model = vcn_model_create();
	vcn_model_get_union(model, model1, model2, 0);
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
	vcn_model_t *model = vcn_model_create();
	vcn_model_get_substraction(model, model1, model2, 0);
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
	vcn_model_t *model = vcn_model_create();
	vcn_model_get_substraction(model, model1, model2, 0);
	TEMPORAL(model); /* TEMPORAL */
	vcn_model_destroy(model1);
	vcn_model_destroy(model2);
	vcn_model_destroy(model);
	return false;
}

static bool check_get_difference_squares(void)
{
	char input_name[256];
	sprintf(input_name, "%s/combining_squares_1.psl", INPUTS_DIR);
	vcn_model_t *model1 = vcn_model_load(input_name);
	sprintf(input_name, "%s/combining_squares_2.psl", INPUTS_DIR);
	vcn_model_t *model2 = vcn_model_load(input_name);
	vcn_model_t *model = vcn_model_create();
	vcn_model_get_difference(model, model1, model2, 0);
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
	vcn_model_t *model = vcn_model_create();
	vcn_model_get_combination(model, model1, model2, 0);
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
	vcn_model_t *model = vcn_model_create();
	vcn_model_get_intersection(model, model1, model2, 0);
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
	vcn_model_t *model = vcn_model_create();
	vcn_model_get_union(model, model1, model2, 0);
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
	vcn_model_t *model = vcn_model_create();
	vcn_model_get_substraction(model, model1, model2, 0);
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
	vcn_model_t *model = vcn_model_create();
	vcn_model_get_substraction(model, model1, model2, 0);
	TEMPORAL(model); /* TEMPORAL */
	vcn_model_destroy(model1);
	vcn_model_destroy(model2);
	vcn_model_destroy(model);
	return false;
}

static bool check_get_difference_cross(void)
{
	char input_name[256];
	sprintf(input_name, "%s/combining_cross_1.psl", INPUTS_DIR);
	vcn_model_t *model1 = vcn_model_load(input_name);
	sprintf(input_name, "%s/combining_cross_2.psl", INPUTS_DIR);
	vcn_model_t *model2 = vcn_model_load(input_name);
	vcn_model_t *model = vcn_model_create();
	vcn_model_get_difference(model, model1, model2, 0);
	TEMPORAL(model); /* TEMPORAL */
	vcn_model_destroy(model1);
	vcn_model_destroy(model2);
	vcn_model_destroy(model);
	return false;
}

static bool check_get_combination_c(void)
{
	char input_name[256];
	sprintf(input_name, "%s/combining_F1.psl", INPUTS_DIR);
	vcn_model_t *model1 = vcn_model_load(input_name);
	sprintf(input_name, "%s/combining_E1.psl", INPUTS_DIR);
	vcn_model_t *model2 = vcn_model_load(input_name);
	vcn_model_t *model = vcn_model_create();
	vcn_model_get_combination(model, model1, model2, 0);
	TEMPORAL(model); /* TEMPORAL */
	vcn_model_destroy(model1);
	vcn_model_destroy(model2);
	vcn_model_destroy(model);
	return false;
}

static bool check_get_intersection_c(void)
{
	char input_name[256];
	sprintf(input_name, "%s/combining_F1.psl", INPUTS_DIR);
	vcn_model_t *model1 = vcn_model_load(input_name);
	sprintf(input_name, "%s/combining_E1.psl", INPUTS_DIR);
	vcn_model_t *model2 = vcn_model_load(input_name);
	vcn_model_t *model = vcn_model_create();
	vcn_model_get_intersection(model, model1, model2, 0);
	TEMPORAL(model); /* TEMPORAL */
	vcn_model_destroy(model1);
	vcn_model_destroy(model2);
	vcn_model_destroy(model);
	return false;
}

static bool check_get_union_c(void)
{
	char input_name[256];
	sprintf(input_name, "%s/combining_F1.psl", INPUTS_DIR);
	vcn_model_t *model1 = vcn_model_load(input_name);
	sprintf(input_name, "%s/combining_E1.psl", INPUTS_DIR);
	vcn_model_t *model2 = vcn_model_load(input_name);
	vcn_model_t *model = vcn_model_create();
	vcn_model_get_union(model, model1, model2, 0);
	TEMPORAL(model); /* TEMPORAL */
	vcn_model_destroy(model1);
	vcn_model_destroy(model2);
	vcn_model_destroy(model);
	return false;
}

static bool check_get_substractionA_c(void)
{
	char input_name[256];
	sprintf(input_name, "%s/combining_F1.psl", INPUTS_DIR);
	vcn_model_t *model1 = vcn_model_load(input_name);
	sprintf(input_name, "%s/combining_E1.psl", INPUTS_DIR);
	vcn_model_t *model2 = vcn_model_load(input_name);
	vcn_model_t *model = vcn_model_create();
	vcn_model_get_substraction(model, model1, model2, 0);
	TEMPORAL(model); /* TEMPORAL */
	vcn_model_destroy(model1);
	vcn_model_destroy(model2);
	vcn_model_destroy(model);
	return false;
}

static bool check_get_substractionB_c(void)
{
	char input_name[256];
	sprintf(input_name, "%s/combining_E1.psl", INPUTS_DIR);
	vcn_model_t *model1 = vcn_model_load(input_name);
	sprintf(input_name, "%s/combining_F1.psl", INPUTS_DIR);
	vcn_model_t *model2 = vcn_model_load(input_name);
	vcn_model_t *model = vcn_model_create();
	vcn_model_get_substraction(model, model1, model2, 0);
	TEMPORAL(model); /* TEMPORAL */
	vcn_model_destroy(model1);
	vcn_model_destroy(model2);
	vcn_model_destroy(model);
	return false;
}

static bool check_get_difference_c(void)
{
	char input_name[256];
	sprintf(input_name, "%s/combining_F1.psl", INPUTS_DIR);
	vcn_model_t *model1 = vcn_model_load(input_name);
	sprintf(input_name, "%s/combining_E1.psl", INPUTS_DIR);
	vcn_model_t *model2 = vcn_model_load(input_name);
	vcn_model_t *model = vcn_model_create();
	vcn_model_get_difference(model, model1, model2, 0);
	TEMPORAL(model); /* TEMPORAL */
	vcn_model_destroy(model1);
	vcn_model_destroy(model2);
	vcn_model_destroy(model);
	return false;
}

static bool check_get_combination_b(void)
{
	char input_name[256];
	sprintf(input_name, "%s/combining_B1.psl", INPUTS_DIR);
	vcn_model_t *model1 = vcn_model_load(input_name);
	sprintf(input_name, "%s/combining_C1.psl", INPUTS_DIR);
	vcn_model_t *model2 = vcn_model_load(input_name);
	vcn_model_t *model = vcn_model_create();
	vcn_model_get_combination(model, model1, model2, 0);
	TEMPORAL(model); /* TEMPORAL */
	vcn_model_destroy(model1);
	vcn_model_destroy(model2);
	vcn_model_destroy(model);
	return false;
}

static bool check_get_intersection_b(void)
{
	char input_name[256];
	sprintf(input_name, "%s/combining_B1.psl", INPUTS_DIR);
	vcn_model_t *model1 = vcn_model_load(input_name);
	sprintf(input_name, "%s/combining_C1.psl", INPUTS_DIR);
	vcn_model_t *model2 = vcn_model_load(input_name);
	vcn_model_t *model = vcn_model_create();
	vcn_model_get_intersection(model, model1, model2, 0);
	TEMPORAL(model); /* TEMPORAL */
	vcn_model_destroy(model1);
	vcn_model_destroy(model2);
	vcn_model_destroy(model);
	return false;
}

static bool check_get_union_b(void)
{
	char input_name[256];
	sprintf(input_name, "%s/combining_B1.psl", INPUTS_DIR);
	vcn_model_t *model1 = vcn_model_load(input_name);
	sprintf(input_name, "%s/combining_C1.psl", INPUTS_DIR);
	vcn_model_t *model2 = vcn_model_load(input_name);
	vcn_model_t *model = vcn_model_create();
	vcn_model_get_union(model, model1, model2, 0);
	TEMPORAL(model); /* TEMPORAL */
	vcn_model_destroy(model1);
	vcn_model_destroy(model2);
	vcn_model_destroy(model);
	return false;
}

static bool check_get_substractionA_b(void)
{
	char input_name[256];
	sprintf(input_name, "%s/combining_B1.psl", INPUTS_DIR);
	vcn_model_t *model1 = vcn_model_load(input_name);
	sprintf(input_name, "%s/combining_C1.psl", INPUTS_DIR);
	vcn_model_t *model2 = vcn_model_load(input_name);
	vcn_model_t *model = vcn_model_create();
	vcn_model_get_substraction(model, model1, model2, 0);
	TEMPORAL(model); /* TEMPORAL */
	vcn_model_destroy(model1);
	vcn_model_destroy(model2);
	vcn_model_destroy(model);
	return false;
}

static bool check_get_substractionB_b(void)
{
	char input_name[256];
	sprintf(input_name, "%s/combining_C1.psl", INPUTS_DIR);
	vcn_model_t *model1 = vcn_model_load(input_name);
	sprintf(input_name, "%s/combining_B1.psl", INPUTS_DIR);
	vcn_model_t *model2 = vcn_model_load(input_name);
	vcn_model_t *model = vcn_model_create();
	vcn_model_get_substraction(model, model1, model2, 0);
	TEMPORAL(model); /* TEMPORAL */
	vcn_model_destroy(model1);
	vcn_model_destroy(model2);
	vcn_model_destroy(model);
	return false;
}

static bool check_get_difference_b(void)
{
	char input_name[256];
	sprintf(input_name, "%s/combining_B1.psl", INPUTS_DIR);
	vcn_model_t *model1 = vcn_model_load(input_name);
	sprintf(input_name, "%s/combining_C1.psl", INPUTS_DIR);
	vcn_model_t *model2 = vcn_model_load(input_name);
	vcn_model_t *model = vcn_model_create();
	vcn_model_get_difference(model, model1, model2, 0);
	TEMPORAL(model); /* TEMPORAL */
	vcn_model_destroy(model1);
	vcn_model_destroy(model2);
	vcn_model_destroy(model);
	return false;
}

static bool check_get_combination_g(void)
{
	char input_name[256];
	sprintf(input_name, "%s/combining_G1.psl", INPUTS_DIR);
	vcn_model_t *model1 = vcn_model_load(input_name);
	sprintf(input_name, "%s/combining_G2.psl", INPUTS_DIR);
	vcn_model_t *model2 = vcn_model_load(input_name);
	vcn_model_t *model = vcn_model_create();
	vcn_model_get_combination(model, model1, model2, 0);
	TEMPORAL(model); /* TEMPORAL */
	vcn_model_destroy(model1);
	vcn_model_destroy(model2);
	vcn_model_destroy(model);
	return false;
}

static bool check_get_intersection_g(void)
{
	char input_name[256];
	sprintf(input_name, "%s/combining_G1.psl", INPUTS_DIR);
	vcn_model_t *model1 = vcn_model_load(input_name);
	sprintf(input_name, "%s/combining_G2.psl", INPUTS_DIR);
	vcn_model_t *model2 = vcn_model_load(input_name);
	vcn_model_t *model = vcn_model_create();
	vcn_model_get_intersection(model, model1, model2, 0);
	TEMPORAL(model); /* TEMPORAL */
	vcn_model_destroy(model1);
	vcn_model_destroy(model2);
	vcn_model_destroy(model);
	return false;
}

static bool check_get_union_g(void)
{
	char input_name[256];
	sprintf(input_name, "%s/combining_G1.psl", INPUTS_DIR);
	vcn_model_t *model1 = vcn_model_load(input_name);
	sprintf(input_name, "%s/combining_G2.psl", INPUTS_DIR);
	vcn_model_t *model2 = vcn_model_load(input_name);
	vcn_model_t *model = vcn_model_create();
	vcn_model_get_union(model, model1, model2, 0);
	TEMPORAL(model); /* TEMPORAL */
	vcn_model_destroy(model1);
	vcn_model_destroy(model2);
	vcn_model_destroy(model);
	return false;
}

static bool check_get_substractionA_g(void)
{
	char input_name[256];
	sprintf(input_name, "%s/combining_G1.psl", INPUTS_DIR);
	vcn_model_t *model1 = vcn_model_load(input_name);
	sprintf(input_name, "%s/combining_G2.psl", INPUTS_DIR);
	vcn_model_t *model2 = vcn_model_load(input_name);
	vcn_model_t *model = vcn_model_create();
	vcn_model_get_substraction(model, model1, model2, 0);
	TEMPORAL(model); /* TEMPORAL */
	vcn_model_destroy(model1);
	vcn_model_destroy(model2);
	vcn_model_destroy(model);
	return false;
}

static bool check_get_substractionB_g(void)
{
	char input_name[256];
	sprintf(input_name, "%s/combining_G2.psl", INPUTS_DIR);
	vcn_model_t *model1 = vcn_model_load(input_name);
	sprintf(input_name, "%s/combining_G1.psl", INPUTS_DIR);
	vcn_model_t *model2 = vcn_model_load(input_name);
	vcn_model_t *model = vcn_model_create();
	vcn_model_get_substraction(model, model1, model2, 0);
	TEMPORAL(model); /* TEMPORAL */
	vcn_model_destroy(model1);
	vcn_model_destroy(model2);
	vcn_model_destroy(model);
	return false;
}

static bool check_get_difference_g(void)
{
	char input_name[256];
	sprintf(input_name, "%s/combining_G1.psl", INPUTS_DIR);
	vcn_model_t *model1 = vcn_model_load(input_name);
	sprintf(input_name, "%s/combining_G2.psl", INPUTS_DIR);
	vcn_model_t *model2 = vcn_model_load(input_name);
	vcn_model_t *model = vcn_model_create();
	vcn_model_get_difference(model, model1, model2, 0);
	TEMPORAL(model); /* TEMPORAL */
	vcn_model_destroy(model1);
	vcn_model_destroy(model2);
	vcn_model_destroy(model);
	return false;
}

static bool check_get_combination_h(void)
{
	char input_name[256];
	sprintf(input_name, "%s/combining_H1.psl", INPUTS_DIR);
	vcn_model_t *model1 = vcn_model_load(input_name);
	sprintf(input_name, "%s/combining_H2.psl", INPUTS_DIR);
	vcn_model_t *model2 = vcn_model_load(input_name);
	vcn_model_t *model = vcn_model_create();
	vcn_model_get_combination(model, model1, model2, 0);
	TEMPORAL(model); /* TEMPORAL */
	vcn_model_destroy(model1);
	vcn_model_destroy(model2);
	vcn_model_destroy(model);
	return false;
}

static bool check_get_intersection_h(void)
{
	char input_name[256];
	sprintf(input_name, "%s/combining_H1.psl", INPUTS_DIR);
	vcn_model_t *model1 = vcn_model_load(input_name);
	sprintf(input_name, "%s/combining_H2.psl", INPUTS_DIR);
	vcn_model_t *model2 = vcn_model_load(input_name);
	vcn_model_t *model = vcn_model_create();
	vcn_model_get_intersection(model, model1, model2, 0);
	TEMPORAL(model); /* TEMPORAL */
	vcn_model_destroy(model1);
	vcn_model_destroy(model2);
	vcn_model_destroy(model);
	return false;
}

static bool check_get_union_h(void)
{
	char input_name[256];
	sprintf(input_name, "%s/combining_H1.psl", INPUTS_DIR);
	vcn_model_t *model1 = vcn_model_load(input_name);
	sprintf(input_name, "%s/combining_H2.psl", INPUTS_DIR);
	vcn_model_t *model2 = vcn_model_load(input_name);
	vcn_model_t *model = vcn_model_create();
	vcn_model_get_union(model, model1, model2, 0);
	TEMPORAL(model); /* TEMPORAL */
	vcn_model_destroy(model1);
	vcn_model_destroy(model2);
	vcn_model_destroy(model);
	return false;
}

static bool check_get_substractionA_h(void)
{
	char input_name[256];
	sprintf(input_name, "%s/combining_H1.psl", INPUTS_DIR);
	vcn_model_t *model1 = vcn_model_load(input_name);
	sprintf(input_name, "%s/combining_H2.psl", INPUTS_DIR);
	vcn_model_t *model2 = vcn_model_load(input_name);
	vcn_model_t *model = vcn_model_create();
	vcn_model_get_substraction(model, model1, model2, 0);
	TEMPORAL(model); /* TEMPORAL */
	vcn_model_destroy(model1);
	vcn_model_destroy(model2);
	vcn_model_destroy(model);
	return false;
}

static bool check_get_substractionB_h(void)
{
	char input_name[256];
	sprintf(input_name, "%s/combining_H2.psl", INPUTS_DIR);
	vcn_model_t *model1 = vcn_model_load(input_name);
	sprintf(input_name, "%s/combining_H1.psl", INPUTS_DIR);
	vcn_model_t *model2 = vcn_model_load(input_name);
	vcn_model_t *model = vcn_model_create();
	vcn_model_get_substraction(model, model1, model2, 0);
	TEMPORAL(model); /* TEMPORAL */
	vcn_model_destroy(model1);
	vcn_model_destroy(model2);
	vcn_model_destroy(model);
	return false;
}

static bool check_get_difference_h(void)
{
	char input_name[256];
	sprintf(input_name, "%s/combining_H1.psl", INPUTS_DIR);
	vcn_model_t *model1 = vcn_model_load(input_name);
	sprintf(input_name, "%s/combining_H2.psl", INPUTS_DIR);
	vcn_model_t *model2 = vcn_model_load(input_name);
	vcn_model_t *model = vcn_model_create();
	vcn_model_get_difference(model, model1, model2, 0);
	TEMPORAL(model); /* TEMPORAL */
	vcn_model_destroy(model1);
	vcn_model_destroy(model2);
	vcn_model_destroy(model);
	return false;
}
