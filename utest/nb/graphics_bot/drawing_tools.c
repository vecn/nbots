#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include <CUnit/Basic.h>

#include "nb/memory_bot.h"
#include "nb/graphics_bot/drawing_tools.h"

static int suite_init(void);
static int suite_clean(void);

static void test_drawing(void);
static void draw_test1(nb_graphics_context_t *g, int w, int h,
		       const void *const data);
static void test_show_text(void);
static void show_text(nb_graphics_context_t *g, int w, int h,
		       const void *const data);
static void test_show_text_tt(void);
static void show_text_tt(nb_graphics_context_t *g, int w, int h,
			 const void *const data);

void cunit_nb_graphics_bot_drawing_tools(void)
{
	CU_pSuite suite =
		CU_add_suite("nb/graphics_bot/drawing_tools.c",
			     suite_init, suite_clean);
	CU_add_test(suite, "drawing()", test_drawing);
	CU_add_test(suite, "show_text()", test_show_text);
	CU_add_test(suite, "show_text_tt()", test_show_text_tt);
}

static int suite_init(void)
{
	return 0;
}

static int suite_clean(void)
{
	return 0;
}


static void test_drawing(void)
{
	/* TEMPORAL FAIL */
	nb_graphics_export("../../../test_drawing.png",
			   200, 180, draw_test1, NULL);
	CU_ASSERT(true);
}

static void draw_test1(nb_graphics_context_t *g, int w, int h,
		       const void *const data)
{
	int scale = 2;
	nb_graphics_move_to(g, scale * 15, scale * 15);
	nb_graphics_line_to(g, scale * 85, scale * 65);

	nb_graphics_move_to(g, scale * 10, scale * 10);
	nb_graphics_line_to(g, scale * 90, scale * 10);
	nb_graphics_line_to(g, scale * 90, scale * 70);
	nb_graphics_line_to(g, scale * 10, scale * 70);
	nb_graphics_close_path(g);
	
	nb_graphics_move_to(g, scale * 80, scale * 20);
	nb_graphics_qcurve_to(g, scale * 20, scale * 60,
			      scale * 80, scale * 60);

	nb_graphics_move_to(g, scale * 75, scale * 25);
	nb_graphics_qrcurve_to(g, scale * 25, scale * 55,
			       scale * 25, scale * 25, 1.0f);

	nb_graphics_move_to(g, scale * 20, scale * 30);
	nb_graphics_curve_to(g, scale * 70, scale * 70,
			     scale * 40, scale * 30,
			     scale * 50, scale * 70);

	nb_graphics_stroke(g);

	nb_graphics_move_to(g, scale * 80, scale * 20);
	nb_graphics_qcurve_to(g, scale * 20, scale * 60,
			      scale * 90, scale * 70);

	nb_graphics_move_to(g, scale * 75, scale * 25);
	nb_graphics_qrcurve_to(g, scale * 25, scale * 55,
			       scale * 25, scale * 25, 3.0f);

	nb_graphics_set_source(g, NB_AQUAMARIN);
	nb_graphics_stroke(g);

	nb_graphics_set_circle(g, scale * 50, scale * 40, scale * 7);
	
	nb_graphics_set_source(g, NB_ROSE);
	nb_graphics_fill_preserve(g);
	nb_graphics_set_source(g, NB_AZURE);
	nb_graphics_stroke(g);

	nb_graphics_set_point(g, scale * 10, scale * 10, scale * 7);
	nb_graphics_set_point(g, scale * 90, scale * 10, scale * 7);
	nb_graphics_set_point(g, scale * 10, scale * 70, scale * 7);
	nb_graphics_set_point(g, scale * 90, scale * 70, scale * 7);
	nb_graphics_set_source(g, NB_ORANGE);
	nb_graphics_stroke(g);

	nb_graphics_palette_t *pal =
		nb_graphics_palette_create_preset(NB_RAINBOW);
	nb_graphics_set_circle(g, scale * 10, scale * 40, scale * 7);	
	nb_graphics_set_source_grad(g, NB_LINEAR,
				    scale * 3, scale * 35,
				    scale * 17, scale * 45,
				    pal);
	nb_graphics_fill(g);
	nb_graphics_palette_destroy(pal);
	
	pal = nb_graphics_palette_create_preset(NB_FRENCH);
	nb_graphics_set_circle(g, scale * 30, scale * 40, scale * 7);	
	nb_graphics_set_source_grad(g, NB_RADIAL,
				    scale * 30, scale * 40,
				    scale * 37, scale * 40,
				    pal);
	nb_graphics_fill(g);
	nb_graphics_palette_destroy(pal);
	
	uint8_t col1[4] = {255, 0, 0, 200};
	uint8_t col2[4] = {0, 255, 0, 255};
	uint8_t col3[4] = {0, 0, 255, 255};
	nb_graphics_set_circle(g, scale * 75, scale * 40, scale * 16);	
	nb_graphics_set_source_trg(g, scale * 75, scale * 28,
				   scale * 64, scale * 45,
				   scale * 86, scale * 45,
				   col1, col2, col3);
	nb_graphics_fill(g);

	nb_graphics_set_source(g, NB_BLACK);

	nb_graphics_move_to(g, scale * 50, scale * 20);
	nb_graphics_line_to(g, scale * 85, scale * 20);
	nb_graphics_set_line_width(g, 1.0);
	nb_graphics_stroke(g);

	nb_graphics_move_to(g, scale * 15, scale * 13);
	nb_graphics_line_to(g, scale * 85, scale * 63);
	nb_graphics_move_to(g, scale * 50, scale * 18);
	nb_graphics_line_to(g, scale * 85, scale * 18);
	nb_graphics_set_line_width(g, 0.8);
	nb_graphics_stroke(g);

	nb_graphics_move_to(g, scale * 15, scale * 11);
	nb_graphics_line_to(g, scale * 85, scale * 61);
	nb_graphics_move_to(g, scale * 50, scale * 16);
	nb_graphics_line_to(g, scale * 85, scale * 16);
	nb_graphics_set_line_width(g, 0.5);
	nb_graphics_stroke(g);

	nb_graphics_move_to(g, scale * 15, scale * 9);
	nb_graphics_line_to(g, scale * 85, scale * 59);
	nb_graphics_move_to(g, scale * 50, scale * 14);
	nb_graphics_line_to(g, scale * 85, scale * 14);
	nb_graphics_set_line_width(g, 0.3);
	nb_graphics_stroke(g);

	nb_graphics_move_to(g, scale * 15, scale * 7);
	nb_graphics_line_to(g, scale * 85, scale * 57);
	nb_graphics_move_to(g, scale * 50, scale * 12);
	nb_graphics_line_to(g, scale * 85, scale * 12);
	nb_graphics_set_line_width(g, 0.1);
	nb_graphics_stroke(g);

	nb_graphics_set_source_rgb(g, 0, 0, 100);

	nb_graphics_move_to(g, scale * 15, scale * 19);
	nb_graphics_line_to(g, scale * 85, scale * 69);
	nb_graphics_move_to(g, scale * 50, scale * 24);
	nb_graphics_line_to(g, scale * 85, scale * 24);
	nb_graphics_set_line_width(g, 1.5);
	nb_graphics_stroke(g);
	
	nb_graphics_move_to(g, scale * 15, scale * 23);
	nb_graphics_line_to(g, scale * 85, scale * 73);
	nb_graphics_move_to(g, scale * 50, scale * 28);
	nb_graphics_line_to(g, scale * 85, scale * 28);
	nb_graphics_set_line_width(g, 2.0);
	nb_graphics_stroke(g);
	
	nb_graphics_move_to(g, scale * 15, scale * 27);
	nb_graphics_line_to(g, scale * 85, scale * 77);
	nb_graphics_move_to(g, scale * 50, scale * 32);
	nb_graphics_line_to(g, scale * 85, scale * 32);
	nb_graphics_set_line_width(g, 2.5);
	nb_graphics_stroke(g);

	nb_graphics_move_to(g, scale * 15, scale * 31);
	nb_graphics_line_to(g, scale * 85, scale * 81);
	nb_graphics_set_line_width(g, 3.0);
	nb_graphics_stroke(g);

	nb_graphics_move_to(g, scale * 15, scale * 35);
	nb_graphics_line_to(g, scale * 85, scale * 85);
	nb_graphics_set_line_width(g, 4.0);
	nb_graphics_stroke(g);

	nb_graphics_move_to(g, scale * 15, scale * 39);
	nb_graphics_line_to(g, scale * 85, scale * 89);
	nb_graphics_set_line_width(g, 5.0);
	nb_graphics_stroke(g);

	nb_graphics_move_to(g, scale * 15, scale * 45);
	nb_graphics_line_to(g, scale * 85, scale * 95);
	nb_graphics_set_line_width(g, 6.0);
	nb_graphics_stroke(g);

	
	nb_graphics_move_to(g, scale * 20, scale * 45);
	nb_graphics_curve_to(g, scale * 70, scale * 85,
			     scale * 40, scale * 45,
			     scale * 50, scale * 85);
	nb_graphics_set_source(g, NB_VIOLET);
	nb_graphics_set_line_width(g, 10.0);
	nb_graphics_stroke(g);
}

static void test_show_text(void)
{
	/* TEMPORAL FAIL */
	nb_graphics_export("../../../test_show_text.png",
			   1000, 500, show_text, NULL);
	CU_ASSERT(true);
}
static void show_text(nb_graphics_context_t *g, int w, int h,
		       const void *const datxa)
{
	const char *label = "Hello, "
		"This test is to show our bitmap font ;)    "
		"1234567890";
	
	nb_graphics_set_font_type(g, "");

	nb_graphics_set_font_size(g, 1);
	nb_graphics_show_text(g, 10, 10, label);

	nb_graphics_set_font_size(g, 2);
	nb_graphics_show_text(g, 10, 20, label);

	nb_graphics_set_font_size(g, 3);
	nb_graphics_show_text(g, 10, 30, label);

	nb_graphics_set_font_size(g, 4);
	nb_graphics_show_text(g, 10, 40, label);

	nb_graphics_set_font_size(g, 5);
	nb_graphics_show_text(g, 10, 50, label);

	nb_graphics_set_font_size(g, 6);
	nb_graphics_show_text(g, 10, 60, label);

	nb_graphics_set_font_size(g, 7);
	nb_graphics_show_text(g, 10, 70, label);

	nb_graphics_set_font_size(g, 8);
	nb_graphics_show_text(g, 10, 80, label);

	nb_graphics_set_font_size(g, 9);
	nb_graphics_show_text(g, 10, 100, label);

	nb_graphics_set_font_size(g, 10);
	nb_graphics_show_text(g, 10, 120, label);

	nb_graphics_set_font_size(g, 11);
	nb_graphics_show_text(g, 10, 140, label);

	nb_graphics_set_font_size(g, 12);
	nb_graphics_show_text(g, 10, 160, label);

	nb_graphics_set_font_size(g, 13);
	nb_graphics_show_text(g, 10, 180, label);

	nb_graphics_set_font_size(g, 14);
	nb_graphics_show_text(g, 10, 200, label);

	nb_graphics_set_font_size(g, 15);
	nb_graphics_show_text(g, 10, 220, label);

	nb_graphics_set_font_size(g, 16);
	nb_graphics_show_text(g, 10, 240, label);

	nb_graphics_set_source(g, NB_BLUE);
	const char *multiline = "Hello world, this is a multiline string\n"
		"to test our library, if anything is wrong please\n"
		"send an email to victorc@cimat.mx ...:\n(;'\"#$%&/=?[]{})";

	nb_graphics_set_font_size(g, 8);
	nb_graphics_show_text(g, 10, 290, multiline);


	nb_graphics_set_font_size(g, 16);
	nb_graphics_text_attr_t attr;
	nb_graphics_get_text_attr(g, multiline, &attr);
	nb_graphics_show_text(g, 500 - attr.width/2,
			      395 + attr.height/2, multiline);
}

static void test_show_text_tt(void)
{
	/* TEMPORAL FAIL */
	nb_graphics_export("../../../test_show_text_tt.png",
			   1000, 500, show_text_tt, NULL);
	CU_ASSERT(true);
}
static void show_text_tt(nb_graphics_context_t *g, int w, int h,
		       const void *const datxa)
{
	const char *label = "Hello world, "
		"This test is to show True type rasterizer "
		"1234567890";
	
	nb_graphics_set_font_type(g, "FreeSerif");

	nb_graphics_set_font_size(g, 1);
	nb_graphics_show_text(g, 10, 10, label);

	nb_graphics_set_font_size(g, 2);
	nb_graphics_show_text(g, 10, 20, label);

	nb_graphics_set_font_size(g, 3);
	nb_graphics_show_text(g, 10, 30, label);

	nb_graphics_set_font_size(g, 4);
	nb_graphics_show_text(g, 10, 40, label);

	nb_graphics_set_font_size(g, 5);
	nb_graphics_show_text(g, 10, 50, label);

	nb_graphics_set_font_size(g, 6);
	nb_graphics_show_text(g, 10, 60, label);

	nb_graphics_set_font_size(g, 7);
	nb_graphics_show_text(g, 10, 70, label);

	nb_graphics_set_font_size(g, 8);
	nb_graphics_show_text(g, 10, 80, label);

	nb_graphics_set_font_size(g, 9);
	nb_graphics_show_text(g, 10, 100, label);

	nb_graphics_set_font_size(g, 10);
	nb_graphics_show_text(g, 10, 120, label);

	nb_graphics_set_font_size(g, 11);
	nb_graphics_show_text(g, 10, 140, label);

	nb_graphics_set_font_size(g, 12);
	nb_graphics_show_text(g, 10, 160, label);

	nb_graphics_set_font_size(g, 13);
	nb_graphics_show_text(g, 10, 180, label);

	nb_graphics_set_font_size(g, 14);
	nb_graphics_show_text(g, 10, 200, label);

	nb_graphics_set_font_size(g, 15);
	nb_graphics_show_text(g, 10, 220, label);

	nb_graphics_set_font_size(g, 16);
	nb_graphics_show_text(g, 10, 240, label);

	nb_graphics_set_source(g, NB_AZURE);
	const char *label_font = "This is a sample "
		"to show a True Type Font 0123456789";

	nb_graphics_set_font_type(g, "FreeMono");
	nb_graphics_set_font_size(g, 20);
	nb_graphics_show_text(g, 10, 270, label_font);

	nb_graphics_set_font_type(g, "FreeSans");
	nb_graphics_set_font_size(g, 20);
	nb_graphics_show_text(g, 10, 300, label_font);

	nb_graphics_set_font_type(g, "FreeSerif");
	nb_graphics_set_font_size(g, 20);
	nb_graphics_show_text(g, 10, 330, label_font);


	nb_graphics_set_source(g, NB_BLUE);
	const char *multiline = "Hello world, this is a multiline string\n"
		"to test our library, if anything is wrong please\n"
		"send an email to victorc@cimat.mx ...:\n(;'\"#$%&/=?[]{})\n"
		"Environment var for TT fonts is NB_FONTS_DIR";

	nb_graphics_set_font_type(g, "FreeSerif");
	nb_graphics_set_font_size(g, 25);
	nb_graphics_text_attr_t attr;
	nb_graphics_get_text_attr(g, multiline, &attr);
	nb_graphics_show_text(g, 500 - attr.width/2,
			      395 + attr.height/2, multiline);
}
