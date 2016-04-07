#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include  "nb/geometric_bot/model/model2D.h"
#include  "nb/geometric_bot/mesh/mesh2D.h"
#include  "nb/geometric_bot/mesh/modules2D/image_density.h"

#include "test_library.h"
#include "test_add.h"

#define INPUTS_DIR "../tests/nb/geometric_bot/mesh/\
modules2D/image_density_UT_inputs"
#define OUTPUT "build"

#include  "nb/geometric_bot/mesh/modules2D/exporter_cairo.h" /* TEMPORAL */
static int TEMPORAL_ = 0; /* TEMPORAL */		      /* TEMPORAL */
static void TEMPORAL(const vcn_mesh_t *const mesh)	      /* TEMPORAL */
{							      /* TEMPORAL */
	char label[100];				      /* TEMPORAL */
	sprintf(label, "%s/TEMP_IMG_%02i.png", OUTPUT,
		TEMPORAL_++);                                 /* TEMPORAL */
	vcn_mesh_save_png(mesh, label, 1000, 800);	      /* TEMPORAL */
}                                                             /* TEMPORAL */

static bool check_set_img_density_jpg_eye(void);
static bool check_set_img_density_jpg_gnome(void);
static bool check_set_img_density_png_jolie(void);
static bool check_set_img_density_jpg_hand(void);
static bool check_set_img_density_jpg_size_const(void);

inline int vcn_test_get_driver_id(void)
{
	return NB_DRIVER_UNIT_TEST;
}

void vcn_test_load_tests(void *tests_ptr)
{
	vcn_test_add(tests_ptr, check_set_img_density_jpg_eye,
		     "Check set_img_density() with a JPG of 'the eye'");
	vcn_test_add(tests_ptr, check_set_img_density_jpg_gnome,
		     "Check set_img_density() with a JPG of 'gnome logo'");
	vcn_test_add(tests_ptr, check_set_img_density_png_jolie,
		     "Check set_img_density() with a PNG of 'Angeline Jolie'");
	vcn_test_add(tests_ptr, check_set_img_density_jpg_hand,
		     "Check set_img_density() with a JPG of 'Hands'");
	vcn_test_add(tests_ptr, check_set_img_density_jpg_size_const,
		     "Check set_img_density() with a JPG with size constraint");
}

static bool check_set_img_density_jpg_eye(void)
{
	char input_name[256];
	sprintf(input_name, "%s/eye.jpg", INPUTS_DIR);
	vcn_image_t *img = vcn_image_create();
	vcn_image_read(img, input_name);	
	vcn_model_t* model = 
		vcn_model_create_rectangle(0.0, 0.0,
					   vcn_image_get_width(img),
					   vcn_image_get_height(img));
	vcn_mesh_t* mesh = vcn_mesh_create();
	vcn_mesh_set_img_density(mesh, img, 0.0);
	vcn_mesh_set_geometric_constraint(mesh,
					  NB_MESH_GEOM_CONSTRAINT_MIN_ANGLE,
					  0.2);
	vcn_mesh_generate_from_model(mesh, model);
	vcn_model_destroy(model);
	vcn_mesh_clear_img_density(mesh);
	vcn_image_destroy(img);
	TEMPORAL(mesh); /* TEMPORAL */
	vcn_mesh_destroy(mesh);
	return false;
}

static bool check_set_img_density_jpg_gnome(void)
{
	char input_name[256];
	sprintf(input_name, "%s/gnome.jpg", INPUTS_DIR);
	vcn_image_t *img = vcn_image_create();
	vcn_image_read(img, input_name);	
	vcn_model_t* model = 
		vcn_model_create_rectangle(0.0, 0.0,
					   vcn_image_get_width(img),
					   vcn_image_get_height(img));
	vcn_mesh_t* mesh = vcn_mesh_create();
	vcn_mesh_set_img_density(mesh, img, 0.0);
	vcn_mesh_generate_from_model(mesh, model);
	vcn_model_destroy(model);
	vcn_mesh_clear_img_density(mesh);
	vcn_image_destroy(img);
	TEMPORAL(mesh); /* TEMPORAL */
	vcn_mesh_destroy(mesh);
	return false;
}

static bool check_set_img_density_png_jolie(void)
{
	char input_name[256];
	sprintf(input_name, "%s/jolie.png", INPUTS_DIR);
	vcn_image_t *img = vcn_image_create();
	vcn_image_read(img, input_name);	
	vcn_model_t* model = 
		vcn_model_create_rectangle(0.0, 0.0,
					   vcn_image_get_width(img),
					   vcn_image_get_height(img));
	vcn_mesh_t* mesh = vcn_mesh_create();
	vcn_mesh_set_img_density(mesh, img, 0.0);
	vcn_mesh_set_geometric_constraint(mesh,
					  NB_MESH_GEOM_CONSTRAINT_MIN_ANGLE,
					  0.2);
	vcn_mesh_generate_from_model(mesh, model);
	vcn_model_destroy(model);
	vcn_mesh_clear_img_density(mesh);
	vcn_image_destroy(img);
	TEMPORAL(mesh); /* TEMPORAL */
	vcn_mesh_destroy(mesh);
	return false;
}

static bool check_set_img_density_jpg_hand(void)
{
	char input_name[256];
	sprintf(input_name, "%s/hand.jpg", INPUTS_DIR);
	vcn_image_t *img = vcn_image_create();
	vcn_image_read(img, input_name);	
	vcn_model_t* model = 
		vcn_model_create_rectangle(0.0, 0.0,
					   vcn_image_get_width(img),
					   vcn_image_get_height(img));
	vcn_mesh_t* mesh = vcn_mesh_create();
	vcn_mesh_set_img_density(mesh, img, 0.0);
	vcn_mesh_set_geometric_constraint(mesh,
					  NB_MESH_GEOM_CONSTRAINT_MIN_ANGLE,
					  0.2);
	vcn_mesh_generate_from_model(mesh, model);
	vcn_model_destroy(model);
	vcn_mesh_clear_img_density(mesh);
	TEMPORAL(mesh); /* TEMPORAL */
	vcn_mesh_destroy(mesh);
	return false;
}

static bool check_set_img_density_jpg_size_const(void)
{
	char input_name[256];
	sprintf(input_name, "%s/color_eye.jpg", INPUTS_DIR);
	vcn_image_t *img = vcn_image_create();
	vcn_image_read(img, input_name);	
	vcn_model_t* model = 
		vcn_model_create_rectangle(0.0, 0.0,
					   vcn_image_get_width(img),
					   vcn_image_get_height(img));
	vcn_mesh_t* mesh = vcn_mesh_create();
	vcn_mesh_set_img_density(mesh, img, 0.0);
	vcn_mesh_set_size_constraint(mesh,
				     NB_MESH_SIZE_CONSTRAINT_MAX_VTX,
				     3000);
	vcn_mesh_set_geometric_constraint(mesh,
					  NB_MESH_GEOM_CONSTRAINT_MIN_ANGLE,
					  0.2);
	vcn_mesh_generate_from_model(mesh, model);
	vcn_model_destroy(model);
	vcn_mesh_clear_img_density(mesh);
	vcn_image_destroy(img);
	TEMPORAL(mesh); /* TEMPORAL */
	vcn_mesh_destroy(mesh);
	return false;
}
