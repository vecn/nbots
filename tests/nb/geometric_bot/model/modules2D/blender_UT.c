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
	vcn_model_export_png(model, label, 1000, 800, false);  /* TEMPORAL */
}                                                             /* TEMPORAL */

static bool check_get_combination(void);

inline int vcn_test_get_driver_id(void)
{
	return NB_DRIVER_UNIT_TEST;
}

void vcn_test_load_tests(void *tests_ptr)
{
	vcn_test_add(tests_ptr, check_get_combination,
		     "Check get_combination()");

}

static bool check_get_combination(void){
	char input_name[256];
	sprintf(input_name, "%s/combining_A1.psl", INPUTS_DIR);
	vcn_model_t *model1 = vcn_model_load(input_name);
	sprintf(input_name, "%s/combining_A2.psl", INPUTS_DIR);
	vcn_model_t *model2 = vcn_model_load(input_name);
	vcn_model_t *model = vcn_model_get_combination(model1, model2, 0);
	TEMPORAL(model); /* TEMPORAL */
	vcn_model_destroy(model1);
	vcn_model_destroy(model2);
	return model;
}
/*
static vcn_model_t* testB(char* name, const char* input_dir){
  sprintf(name, "Combining test B");
  char input_name[256];
  sprintf(input_name, "%s/combining_A1.psl", input_dir);
  vcn_model_t *model1 = vcn_model_load(input_name);
  sprintf(input_name, "%s/combining_B1.psl", input_dir);
  vcn_model_t *model2 = vcn_model_load(input_name);
  vcn_model_t *model = function_to_test(model1, model2, 0);
  vcn_model_destroy(model1);
  vcn_model_destroy(model2);
  return model;
}

static vcn_model_t* testC(char* name, const char* input_dir){
  sprintf(name, "Combining test C");
  char input_name[256];
  sprintf(input_name, "%s/combining_A1.psl", input_dir);
  vcn_model_t *model1 = vcn_model_load(input_name);
  sprintf(input_name, "%s/combining_C1.psl", input_dir);
  vcn_model_t *model2 = vcn_model_load(input_name);
  vcn_model_t *model = function_to_test(model1, model2, 0);
  vcn_model_destroy(model1);
  vcn_model_destroy(model2);
  return model;
}

static vcn_model_t* testD(char* name, const char* input_dir){
  sprintf(name, "Combining test D");
  char input_name[256];
  sprintf(input_name, "%s/combining_D1.psl", input_dir);
  vcn_model_t *model1 = vcn_model_load(input_name);
  sprintf(input_name, "%s/combining_D2.psl", input_dir);
  vcn_model_t *model2 = vcn_model_load(input_name);
  vcn_model_t *model = function_to_test(model1, model2, 0);
  vcn_model_destroy(model1);
  vcn_model_destroy(model2);
  return model;
}

static vcn_model_t* testE(char* name, const char* input_dir){
  sprintf(name, "Combining test E");
  char input_name[256];
  sprintf(input_name, "%s/combining_A1.psl", input_dir);
  vcn_model_t *model1 = vcn_model_load(input_name);
  sprintf(input_name, "%s/combining_E1.psl", input_dir);
  vcn_model_t *model2 = vcn_model_load(input_name);
  vcn_model_t *model = function_to_test(model1, model2, 0);
  vcn_model_destroy(model1);
  vcn_model_destroy(model2);
  return model;
}

static vcn_model_t* testF(char* name, const char* input_dir){
  sprintf(name, "Combining test F");
  char input_name[256];
  sprintf(input_name, "%s/combining_A1.psl", input_dir);
  vcn_model_t *model1 = vcn_model_load(input_name);
  sprintf(input_name, "%s/combining_F1.psl", input_dir);
  vcn_model_t *model2 = vcn_model_load(input_name);
  vcn_model_t *model = function_to_test(model1, model2, 0);
  vcn_model_destroy(model1);
  vcn_model_destroy(model2);
  return model;
}

static vcn_model_t* testG(char* name, const char* input_dir){
  sprintf(name, "Combining test G");
  char input_name[256];
  sprintf(input_name, "%s/combining_G1.psl", input_dir);
  vcn_model_t *model1 = vcn_model_load(input_name);
  sprintf(input_name, "%s/combining_G2.psl", input_dir);
  vcn_model_t *model2 = vcn_model_load(input_name);
  vcn_model_t *model = function_to_test(model1, model2, 0);
  vcn_model_destroy(model1);
  vcn_model_destroy(model2);
  return model;
}

static vcn_model_t* testH(char* name, const char* input_dir){
  sprintf(name, "Combining test H");
  char input_name[256];
  sprintf(input_name, "%s/combining_H1.psl", input_dir);
  vcn_model_t *model1 = vcn_model_load(input_name);
  sprintf(input_name, "%s/combining_H2.psl", input_dir);
  vcn_model_t *model2 = vcn_model_load(input_name);
  vcn_model_t *model = function_to_test(model1, model2, 0);
  vcn_model_destroy(model1);
  vcn_model_destroy(model2);
  return model;
}
*/
