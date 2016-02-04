#include <stdbool.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include  "vcn/geometric_bot/mesh/modules2D/image_density.h"

#include "test_library.h"
#include "test_add.h"

#define INPUTS_DIR "../tests/vcn/geometric_bot/mesh/\
modules2D/image_density_UT_inputs"

static bool check_some(void);

inline int vcn_test_get_driver_id(void)
{
	return VCN_DRIVER_UNIT_TEST;
}

void vcn_test_load_tests(void *tests_ptr)
{
	vcn_test_add(tests_ptr, check_some,
		     "Check some()");
}

/**
 * @brief Test the PNG format in a color picture.
 */
/*
static vcn_mesh_t* test17(char* name, const char* input_dir)
{
  char input_name[256];
  sprintf(input_name, "%s/eye_raw.jpg", input_dir);
  sprintf(name, "Eye picture");
  vcn_density_img_t* data = vcn_density_img_create(input_name, 1.0, 
						   0.0, 0.0, 0.2);
  vcn_model_t* model = 
    vcn_model_create_rectangle(0.0, 0.0,
			       vcn_density_img_get_width(data),
			       vcn_density_img_get_height(data));
  vcn_mesh_t* mesh = 
    vcn_mesh_create_from_model(model, 0, 0, 0.2, VCN_DENSITY_IMG, data);
  vcn_model_destroy(model);
  vcn_density_img_destroy(data);

  return mesh;
}
*/

/**
 * @brief Test the JPEG format.
 */
/*
static vcn_mesh_t* test18(char* name, const char* input_dir)
{
  char input_name[256];
  sprintf(input_name, "%s/gnome.jpg", input_dir);
  sprintf(name, "Gnome logo");
  vcn_density_img_t* data = vcn_density_img_create(input_name, 1.0, 
						   0.0, 0.0, 1.0);
  vcn_model_t* model = 
    vcn_model_create_rectangle(0.0, 0.0,
			       vcn_density_img_get_width(data),
			       vcn_density_img_get_height(data));
  vcn_mesh_t* mesh = 
    vcn_mesh_create_from_model(model, 0, 0,
			       VCN_ANGLE_MAX,
			       VCN_DENSITY_IMG, data);
  vcn_model_destroy(model);
  vcn_density_img_destroy(data);

  return mesh;
}
*/

/**
 * @brief Test a grayscale PNG image.
 */
/*
static vcn_mesh_t* test19(char* name, const char* input_dir)
{
  char input_name[256];
  sprintf(input_name, "%s/women.png", input_dir);
  sprintf(name, "Women picture");

  vcn_density_img_t* data = vcn_density_img_create(input_name, 1.0, 
						   0.0, 0.0, 2.0);
  vcn_model_t* model = 
    vcn_model_create_rectangle(0.0, 0.0,
			       vcn_density_img_get_width(data),
			       vcn_density_img_get_height(data));
  vcn_mesh_t* mesh = 
    vcn_mesh_create_from_model(model, 0, 0, 0.2, VCN_DENSITY_IMG, data);
  vcn_model_destroy(model);
  vcn_density_img_destroy(data);

  return mesh;
}
*/

/**
 * @brief Test JPG black and white image.
 */
/*
static vcn_mesh_t* test20(char* name, const char* input_dir)
{
  char input_name[256];
  sprintf(input_name, "%s/hand.jpg", input_dir);
  sprintf(name, "Hand and baby hand");

  vcn_density_img_t* data = vcn_density_img_create(input_name, 1.0, 
						   0.0, 0.0, 0.2);
  vcn_model_t* model =
    vcn_model_create_rectangle(0.0, 0.0,
			       vcn_density_img_get_width(data),
			       vcn_density_img_get_height(data));
  vcn_mesh_t* mesh = 
    vcn_mesh_create_from_model(model, 0, 0, 0.2, VCN_DENSITY_IMG, data);
  vcn_model_destroy(model);
  vcn_density_img_destroy(data);

  return mesh;
}
*/
