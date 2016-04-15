#include <stdbool.h>
#include <stdlib.h>
#include <alloca.h>

#include <CUnit/Basic.h>

#include "nb/geometric_bot/point2D.h"

#define OUTPUT_DIR "../../../"

#include "nb/geometric_bot/mesh/modules2D/exporter_cairo.h"   /* TEMPORAL */
static int TEMPORAL_ = 0; /* TEMPORAL */		      /* TEMPORAL */
#include <stdio.h>                                           /* TEMPORAL */
static void TEMPORAL(const nb_mshquad_t *const quad)	      /* TEMPORAL */
{							      /* TEMPORAL */
	char label[100];				      /* TEMPORAL */
	sprintf(label, "%s/QUAD_%02i.png", OUTPUT_DIR,
		TEMPORAL_++);	                              /* TEMPORAL */
	nb_mshquad_export_png(quad, label, 1000, 800);	      /* TEMPORAL */
}                                                             /* TEMPORAL */

static int suite_init(void);
static int suite_clean(void);

static void test_load_from_mesh(void);

void cunit_nb_geometric_bot_mesh_elements2D_quad(void)
{
	CU_pSuite suite = CU_add_suite("nb/geometric_bot/mesh/elements2D/quad.c",
				       suite_init, suite_clean);
	CU_add_test(suite, "load_from_mesh()", test_load_from_mesh);
}

static int suite_init(void)
{
	return 0;
}

static int suite_clean(void)
{
	return 0;
}

static void test_load_from_mesh(void)
{
	vcn_model_t *model = vcn_model_create_polygon(20, 0, 0, 100);
	vcn_mesh_t* mesh = vcn_mesh_create();
	vcn_mesh_set_geometric_constraint(mesh,
					 NB_MESH_GEOM_CONSTRAINT_MAX_EDGE_LENGTH,
					 1.0);
	vcn_mesh_generate_from_model(mesh, model);
	vcn_model_destroy(model);
	
	uint32_t size = nb_mshquad_get_memsize();
	nb_mshquad_t *quad = alloca(size);
	nb_mshquad_init(quad);
	nb_mshquad_load_from_mesh(quad, mesh);

	vcn_mesh_destroy(mesh);
	TEMPORAL(quad); /* TEMPORAL */

	nb_mshquad_finish(quad);
	CU_ASSERT(false);
}
