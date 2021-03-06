#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "cunit/Basic.h"

#include "nb/memory_bot.h"
#include "nb/geometric_bot.h"

#define INPUTS_DIR "../utest/sources/nb/geometric_bot/mesh/modules2D/image_density_inputs"

static int suite_init(void);
static int suite_clean(void);

static void test_set_img_density_jpg_eye(void);
static void test_set_img_density_jpg_gnome(void);
static void test_set_img_density_png_jolie(void);
static void test_set_img_density_jpg_hand(void);
static void test_set_img_density_jpg_size_const(void);

void cunit_nb_geometric_bot_tessellator2D_image_density(void)
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

	uint32_t memsize = nb_image_get_memsize();
	nb_image_t *img = nb_soft_allocate_mem(memsize);
	
	nb_image_init(img);
	nb_image_read(img, input_name);	
	nb_model_t* model = 
		nb_model_create_rectangle(0.0, 0.0,
					   nb_image_get_width(img),
					   nb_image_get_height(img));
	nb_tessellator2D_t* mesh = nb_tessellator2D_create();
	nb_tessellator2D_set_img_density(mesh, img, 0.0);
	nb_tessellator2D_set_geometric_constraint(mesh,
					  NB_MESH_GEOM_CONSTRAINT_MIN_ANGLE,
					  0.2);
	nb_tessellator2D_generate_from_model(mesh, model);
	nb_model_destroy(model);
	nb_tessellator2D_clear_img_density(mesh);
	nb_image_finish(img);

	uint32_t N_trg = nb_tessellator2D_get_N_trg(mesh);
	uint32_t N_edge = nb_tessellator2D_get_N_edg(mesh);
	nb_tessellator2D_destroy(mesh);
	/* TEMPORAL FAIL: Produce different triangles each time */
	CU_ASSERT(3600 < N_trg && 3900 > N_trg);
	CU_ASSERT(5500 < N_edge && 5800 > N_edge);

	nb_soft_free_mem(memsize, img);
}

static void test_set_img_density_jpg_gnome(void)
{
	char input_name[256];
	sprintf(input_name, "%s/gnome.jpg", INPUTS_DIR);

	uint32_t memsize = nb_image_get_memsize();
	nb_image_t *img = nb_soft_allocate_mem(memsize);
	
	nb_image_init(img);
	nb_image_read(img, input_name);	
	nb_model_t* model = 
		nb_model_create_rectangle(0.0, 0.0,
					   nb_image_get_width(img),
					   nb_image_get_height(img));
	nb_tessellator2D_t* mesh = nb_tessellator2D_create();
	nb_tessellator2D_set_img_density(mesh, img, 0.0);
	nb_tessellator2D_generate_from_model(mesh, model);
	nb_model_destroy(model);
	nb_tessellator2D_clear_img_density(mesh);
	nb_image_finish(img);
	uint32_t N_trg = nb_tessellator2D_get_N_trg(mesh);
	uint32_t N_edge = nb_tessellator2D_get_N_edg(mesh);
	nb_tessellator2D_destroy(mesh);
	/* TEMPORAL FAIL: Produce different triangles each time */
	CU_ASSERT(12000 < N_trg && 12900 > N_trg);
	CU_ASSERT(18000 < N_edge && 19500 > N_edge);

	nb_soft_free_mem(memsize, img);
}

static void test_set_img_density_png_jolie(void)
{
	char input_name[256];
	sprintf(input_name, "%s/jolie.png", INPUTS_DIR);

	uint32_t memsize = nb_image_get_memsize();
	nb_image_t *img = nb_soft_allocate_mem(memsize);
	
	nb_image_init(img);
	nb_image_read(img, input_name);	
	nb_model_t* model = 
		nb_model_create_rectangle(0.0, 0.0,
					   nb_image_get_width(img),
					   nb_image_get_height(img));
	nb_tessellator2D_t* mesh = nb_tessellator2D_create();
	nb_tessellator2D_set_img_density(mesh, img, 0.0);
	nb_tessellator2D_set_geometric_constraint(mesh,
					  NB_MESH_GEOM_CONSTRAINT_MIN_ANGLE,
					  0.2);
	nb_tessellator2D_generate_from_model(mesh, model);
	nb_model_destroy(model);
	nb_tessellator2D_clear_img_density(mesh);
	nb_image_finish(img);
	uint32_t N_trg = nb_tessellator2D_get_N_trg(mesh);
	uint32_t N_edge = nb_tessellator2D_get_N_edg(mesh);
	nb_tessellator2D_destroy(mesh);
	/* TEMPORAL FAIL: Produce different triangles each time */
	CU_ASSERT(11300 < N_trg && 11900 > N_trg);
	CU_ASSERT(17100 < N_edge && 17900 > N_edge);

	nb_soft_free_mem(memsize, img);
}

static void test_set_img_density_jpg_hand(void)
{
	char input_name[256];
	sprintf(input_name, "%s/hand.jpg", INPUTS_DIR);

	uint32_t memsize = nb_image_get_memsize();
	nb_image_t *img = nb_soft_allocate_mem(memsize);
	
	nb_image_init(img);
	nb_image_read(img, input_name);	
	nb_model_t* model = 
		nb_model_create_rectangle(0.0, 0.0,
					   nb_image_get_width(img),
					   nb_image_get_height(img));
	nb_tessellator2D_t* mesh = nb_tessellator2D_create();
	nb_tessellator2D_set_img_density(mesh, img, 0.0);
	nb_tessellator2D_set_geometric_constraint(mesh,
					  NB_MESH_GEOM_CONSTRAINT_MIN_ANGLE,
					  0.2);
	nb_tessellator2D_generate_from_model(mesh, model);
	nb_model_destroy(model);
	nb_tessellator2D_clear_img_density(mesh);
	nb_image_finish(img);
	uint32_t N_trg = nb_tessellator2D_get_N_trg(mesh);
	uint32_t N_edge = nb_tessellator2D_get_N_edg(mesh);
	nb_tessellator2D_destroy(mesh);
	/* TEMPORAL FAIL: Produce different triangles each time */
	CU_ASSERT(13200 < N_trg && 13700 > N_trg);
	CU_ASSERT(20000 < N_edge && 21000 > N_edge);

	nb_soft_free_mem(memsize, img);
}

static void test_set_img_density_jpg_size_const(void)
{
	char input_name[256];
	sprintf(input_name, "%s/color_eye.jpg", INPUTS_DIR);

	uint32_t memsize = nb_image_get_memsize();
	nb_image_t *img = nb_soft_allocate_mem(memsize);
	
	nb_image_init(img);
	nb_image_read(img, input_name);	
	nb_model_t* model = 
		nb_model_create_rectangle(0.0, 0.0,
					   nb_image_get_width(img),
					   nb_image_get_height(img));
	nb_tessellator2D_t* mesh = nb_tessellator2D_create();
	nb_tessellator2D_set_img_density(mesh, img, 0.0);
	nb_tessellator2D_set_size_constraint(mesh,
				     NB_MESH_SIZE_CONSTRAINT_MAX_VTX,
				     3000);
	nb_tessellator2D_set_geometric_constraint(mesh,
					  NB_MESH_GEOM_CONSTRAINT_MIN_ANGLE,
					  0.2);
	nb_tessellator2D_generate_from_model(mesh, model);
	nb_model_destroy(model);
	nb_tessellator2D_clear_img_density(mesh);
	nb_image_finish(img);
	uint32_t N_trg = nb_tessellator2D_get_N_trg(mesh);
	uint32_t N_edge = nb_tessellator2D_get_N_edg(mesh);
	nb_tessellator2D_destroy(mesh);
	/* TEMPORAL FAIL: Produce different triangles each time */
	CU_ASSERT(5800 < N_trg && 6000 > N_trg);
	CU_ASSERT(8800 < N_edge && 10001 > N_edge);

	nb_soft_free_mem(memsize, img);
}
