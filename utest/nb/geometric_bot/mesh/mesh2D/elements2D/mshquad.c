#include <stdbool.h>
#include <stdlib.h>

#include <CUnit/Basic.h>

#include "nb/memory_bot.h"
#include "nb/geometric_bot.h"

static int suite_init(void);
static int suite_clean(void);

static void test_load_from_tessellator2D(void);

void cunit_nb_mshquad(void)
{
	CU_pSuite suite = CU_add_suite("nb/geometric_bot/mesh/"\
				       "mesh2D/elements2D/mshquad.c",
				       suite_init, suite_clean);
	CU_add_test(suite, "load_from_mesh()", test_load_from_tessellator2D);
}

static int suite_init(void)
{
	return 0;
}

static int suite_clean(void)
{
	return 0;
}

static void test_load_from_tessellator2D(void)
{
	nb_model_t *model = nb_model_create_polygon(1, 0, 0, 6);
	nb_tessellator2D_t* mesh = nb_tessellator2D_create();
	nb_tessellator2D_set_geometric_constraint(mesh,
					 NB_MESH_GEOM_CONSTRAINT_MAX_EDGE_LENGTH,
					 0.1);
	nb_tessellator2D_generate_from_model(mesh, model);
	nb_model_destroy(model);
	
	uint32_t size = nb_mshquad_get_memsize();
	nb_mshquad_t *quad = nb_allocate_on_stack(size);
	nb_mshquad_init(quad);
	nb_mshquad_load_from_tessellator2D(quad, mesh);
	nb_tessellator2D_destroy(mesh);
	
	/* TEMPORAL FAIL: Produce different triangles each time */
	uint32_t N_elems = nb_mshquad_get_N_elems(quad);
	uint32_t N_edg = nb_mshquad_get_N_edges(quad);
	CU_ASSERT(600 < N_elems && 760 > N_elems);
	CU_ASSERT(1200 < N_edg && 1400 > N_edg);

	nb_mshquad_finish(quad);
}
