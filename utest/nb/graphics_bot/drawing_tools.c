#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <alloca.h>

#include <CUnit/Basic.h>

#include "nb/graphics_bot/drawing_tools.h"

static int suite_init(void);
static int suite_clean(void);

static void draw_test1(nb_graphics_context_t *g, int w, int h,
		       const void *const data);
static void test_drawing(void);

void cunit_nb_graphics_bot_drawing_tools(void)
{
	CU_pSuite suite =
		CU_add_suite("nb/graphics_bot/drawing_tools.c", suite_init, suite_clean);
	CU_add_test(suite, "drawing()", test_drawing);
}

static int suite_init(void)
{
	return 0;
}

static int suite_clean(void)
{
	return 0;
}

static void draw_test1(nb_graphics_context_t *g, int w, int h,
		       const void *const data)
{
	nb_graphics_move_to(g, 15, 15);
	nb_graphics_line_to(g, 85, 65);

	nb_graphics_move_to(g, 10, 10);
	nb_graphics_line_to(g, 90, 10);
	nb_graphics_line_to(g, 90, 70);
	nb_graphics_line_to(g, 10, 70);
	nb_graphics_close_path(g);
	
	nb_graphics_move_to(g, 80, 20);
	nb_graphics_qcurve_to(g, 20, 60, 80, 60);

	nb_graphics_move_to(g, 75, 25);
	nb_graphics_qrcurve_to(g, 25, 55, 25, 25, 1.0f);

	nb_graphics_move_to(g, 20, 30);
	nb_graphics_curve_to(g, 70, 70, 40, 30, 50, 70);

	nb_graphics_stroke(g);

	nb_graphics_move_to(g, 80, 20);
	nb_graphics_qcurve_to(g, 20, 60, 90, 70);

	nb_graphics_move_to(g, 75, 25);
	nb_graphics_qrcurve_to(g, 25, 55, 25, 25, 3.0f);

	nb_graphics_set_source(g, NB_AQUAMARIN);
	nb_graphics_stroke(g);

	nb_graphics_set_circle(g, 50, 40, 7);
	
	nb_graphics_set_source(g, NB_ROSE);
	nb_graphics_fill_preserve(g);
	nb_graphics_set_source(g, NB_AZURE);
	nb_graphics_stroke(g);

	nb_graphics_set_point(g, 10, 10, 7);
	nb_graphics_set_point(g, 90, 10, 7);
	nb_graphics_set_point(g, 10, 70, 7);
	nb_graphics_set_point(g, 90, 70, 7);
	nb_graphics_set_source(g, NB_ORANGE);
	nb_graphics_stroke(g);

	nb_graphics_palette_t *pal =
		nb_graphics_palette_create_preset(NB_RAINBOW);
	nb_graphics_set_circle(g, 10, 40, 7);	
	nb_graphics_set_source_grad(g, NB_LINEAR,
				    3, 40, 17, 40,
				    pal);
	nb_graphics_fill(g);
	nb_graphics_palette_destroy(pal);
	
	pal = nb_graphics_palette_create_preset(NB_FRENCH);
	nb_graphics_set_circle(g, 30, 40, 7);	
	nb_graphics_set_source_grad(g, NB_RADIAL,
				    30, 40, 37, 40,
				    pal);
	nb_graphics_fill(g);
	nb_graphics_palette_destroy(pal);
	
	uint8_t col1[4] = {255, 0, 0, 200};
	uint8_t col2[4] = {0, 255, 0, 255};
	uint8_t col3[4] = {0, 0, 255, 255};
	nb_graphics_set_circle(g, 75, 40, 16);	
	nb_graphics_set_source_trg(g,
				   75, 28, 64, 45, 86, 45,
				   col1, col2, col3);
	nb_graphics_fill(g);
}

static void test_drawing(void)
{
	nb_graphics_export("../../../test_drawing.png", 100, 80, draw_test1,
			   NULL);
	CU_ASSERT(true);
}
