#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include <CUnit/Basic.h>

#include "nb/geometric_bot.h"
#include "nb/graphics_bot.h"

#define INPUTS_DIR "../../../../utest/nb/geometric_bot/mesh/modules2D/image_density_inputs"
#define OUTPUT "../../../"

static int TEMPORAL_ = 0; /* TEMPORAL */		      /* TEMPORAL */
static void TEMPORAL(const vcn_mesh_t *const mesh)	      /* TEMPORAL */
{							      /* TEMPORAL */
	char label[100];				      /* TEMPORAL */
	sprintf(label, "%s/TEMP_IMG_%02i.png", OUTPUT,
		TEMPORAL_++);                                 /* TEMPORAL */
	vcn_mesh_save_png(mesh, label, 1000, 800);	      /* TEMPORAL */
}                                                             /* TEMPORAL */

static int suite_init(void);
static int suite_clean(void);

static void test_set_img_density_jpg_eye(void);
static void test_set_img_density_jpg_gnome(void);
static void test_set_img_density_png_jolie(void);
static void test_set_img_density_jpg_hand(void);
static void test_set_img_density_jpg_size_const(void);

void cunit_nb_geometric_bot_mesh2D_image_density(void)
{
	CU_pSuite suite =
		CU_add_suite("nb/geometric_bot/mesh/modules2D/image_density.c",
			     suite_init, suite_clean);
	CU_add_test(suite, "set_img_density() with a JPG of 'the eye'",
		    test_set_img_density_jpg_eye);
	CU_add_test(suite, "set_img_density() with a JPG of 'gnome logo'",
		    test_set_img_density_jpg_gnome);
	CU_add_test(suite, "set_img_density() with a PNG of 'Angeline Jolie'",
		    test_set_img_density_png_jolie);
	CU_add_test(suite, "set_img_density() with a JPG of 'Hands'",
		    test_set_img_density_jpg_hand);
	CU_add_test(suite,
		    "set_img_density() with a JPG with size constraint",
		    test_set_img_density_jpg_size_const);
}

static int suite_init(void)
{
	return 0;
}

static int suite_clean(void)
{
	return 0;
}

static void test_set_img_density_jpg_eye(void)
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
	CU_ASSERT(false);
}

static void test_set_img_density_jpg_gnome(void)
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
	CU_ASSERT(false);
}

static void test_set_img_density_png_jolie(void)
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
	CU_ASSERT(false);
}

static void test_set_img_density_jpg_hand(void)
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
	CU_ASSERT(false);
}

static void test_set_img_density_jpg_size_const(void)
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
	CU_ASSERT(false);
}
