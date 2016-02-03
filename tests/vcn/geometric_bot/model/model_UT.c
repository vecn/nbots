#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "vcn/generic_dst.h"
#include "vcn/ugrid_dst.h"
#include "vcn/sparse_bot.h"
#include "vcn/geometric_bot.h"
#include "vcn/geometric_bot-cairo.h"

static vcn_model_t* testA (char* name, const char* input_dir);
static vcn_model_t* testB (char* name, const char* input_dir);
static vcn_model_t* testC (char* name, const char* input_dir);
static vcn_model_t* testD (char* name, const char* input_dir);
static vcn_model_t* testE (char* name, const char* input_dir);
static vcn_model_t* testF (char* name, const char* input_dir);
static vcn_model_t* testG (char* name, const char* input_dir);
static vcn_model_t* testG (char* name, const char* input_dir);
static vcn_model_t* testH (char* name, const char* input_dir);

static vcn_model_t* (*function_to_test)(const vcn_model_t *const model1,
					const vcn_model_t *const model2,
					double min_length_x_segment);

int main(int argc, char* argv[]){
  if(argc < 3){
    printf("The program must receive:\n");
    printf("   > Directory with the inputs\n");
    printf("   > Directory to save the outputs\n");
    printf("   > [OPTIONAL: Id of operation]\n");
    printf("   > [OPTIONAL: Id of the test to perform]\n\n");
    return 0;
  }

  uint operation_id = 0;
  if(argc > 3)
    operation_id = atoi(argv[3]);

  uint test_id = 0;
  if(argc > 4)
    test_id = atoi(argv[4]);
  
  if(operation_id == 1){
    function_to_test = vcn_model_get_intersection;
    printf("Operation to test: INTERSECTION\n");
  }else if(operation_id == 2){
    function_to_test = vcn_model_get_union;
    printf("Operation to test: UNION\n");
  }else{
    function_to_test = vcn_model_get_combination;
    printf("Operation to test: COMBINATION\n");
  }
  
  const uint N_tests = 8;
  vcn_model_t* (*test[8])(char*, const char*);
  test[0] = testA;
  test[1] = testB;
  test[2] = testC;
  test[3] = testD;
  test[4] = testE;
  test[5] = testF;
  test[6] = testG;
  test[7] = testH;

  /* Execute tests */
  char name[100];

  for(register uint i = 0; i < N_tests; i++){
    if(i != test_id - 1 && test_id != 0) 
      continue;

    printf("Computing next test...");
    fflush(stdout);

    /* Run test */
    vcn_model_t*  model = test[i](name, argv[1]);

    if(model == NULL){
      printf("\r%i: %s is NULL             \n", i+1, name);
      continue;
    }

    int error = vcn_model_verify_consistence(model, NULL);
    if(error != 0){
      printf("\r%i: %s is not coherent (%i).             \n", 
	     i+1, name, error);
      continue;
    }

    /* Show details */
    printf("\r%i: %s              \n", i+1, name);
    printf("    Vertices: %i\n", vcn_model_get_number_of_vertices(model));
    printf("    Segments: %i\n", vcn_model_get_number_of_edges(model));
    printf("       Holes: %i\n", vcn_model_get_number_of_hole_seeds(model));

    /* Visualize Mesh */
    vcn_mesh_t* mesh = 
      vcn_mesh_create_from_model(model, vcn_model_get_number_of_vertices(model),
				 0, 0.0, NULL, NULL);
    vcn_msh3trg_t* msh3trg = 
      vcn_mesh_get_msh3trg(mesh, false, true, false, true, true, NULL);

    sprintf(name, "%s/test%i_mesh.png", argv[2], i+1);
    vcn_msh3trg_save_png(msh3trg, name, 1000, 800,
			 NULL, NULL, NULL, NULL, 0.5, 1.0);

    /* Free memory */
    vcn_msh3trg_destroy(msh3trg);
    vcn_mesh_destroy(mesh);
    vcn_model_destroy(model);
  }

  /* Successful finish */
  return 0;
}

static vcn_model_t* testA(char* name, const char* input_dir){
  sprintf(name, "Combining test A");
  char input_name[256];
  sprintf(input_name, "%s/combining_A1.psl", input_dir);
  vcn_model_t *model1 = vcn_model_load(input_name);
  sprintf(input_name, "%s/combining_A2.psl", input_dir);
  vcn_model_t *model2 = vcn_model_load(input_name);
  vcn_model_t *model = function_to_test(model1, model2, 0);
  vcn_model_destroy(model1);
  vcn_model_destroy(model2);
  return model;
}

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
