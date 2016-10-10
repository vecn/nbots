#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include <CUnit/Basic.h>

#include "nb/memory_bot.h"
#include "nb/image_bot/image.h"

#define INPUTS_DIR "../../../../utest/nb/image_bot/image_inputs"


static int suite_init(void);
static int suite_clean(void);

static void test_render_ascii(void);

void cunit_nb_image_bot_image(void)
{
	CU_pSuite suite =
		CU_add_suite("nb/image_bot/image.c", suite_init, suite_clean);
	CU_add_test(suite, "render_ascii()", test_render_ascii);
}

static int suite_init(void)
{
	return 0;
}

static int suite_clean(void)
{
	return 0;
}

static void test_render_ascii(void)
{
	char filename[100];
	sprintf(filename, "%s/%s", INPUTS_DIR, "jessica.jpg");
	
	nb_image_t *img = nb_allocate_on_stack(nb_image_get_memsize());
	nb_image_init(img);
	nb_image_read(img, filename);
	nb_image_write_ascii(img, "../../../jessica030.txt", 30);
	nb_image_write_ascii(img, "../../../jessica050.txt", 50);
	nb_image_write_ascii(img, "../../../jessica075.txt", 75);
	nb_image_write_ascii(img, "../../../jessica100.txt", 100);
	nb_image_finish(img);
	/* TEMPORAL: this test is not checked */
	CU_ASSERT(true);
}
