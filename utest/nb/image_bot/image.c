#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <alloca.h>

#include <CUnit/Basic.h>

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
	sprintf(filename, "%s/%s", INPUTS_DIR, "ernesto.jpg");
	
	vcn_image_t *img = alloca(vcn_image_get_memsize());
	vcn_image_init(img);
	vcn_image_read(img, filename);
	vcn_image_write_ascii(img, "../../../ernesto030.txt", 30);
	vcn_image_write_ascii(img, "../../../ernesto050.txt", 50);
	vcn_image_write_ascii(img, "../../../ernesto075.txt", 75);
	vcn_image_write_ascii(img, "../../../ernesto100.txt", 100);
	vcn_image_write_ascii(img, "../../../ernesto150.txt", 150);
	vcn_image_finish(img);
	CU_ASSERT(true);
}
