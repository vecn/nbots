#include <stdbool.h>
#include <stdlib.h>
#include <alloca.h>

#include <CUnit/Basic.h>

#include "nb/graphics_bot.h"
#include "nb/geometric_bot.h"

static int suite_init(void);
static int suite_clean(void);

static void test_load_from_mesh(void);

void cunit_nb_geometric_bot_mesh_elements2D_triangles(void)
{
	CU_pSuite suite = CU_add_suite("nb/geometric_bot/mesh/elements2D/triangles.c",
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
	nb_model_t *model  = alloca(vcn_model_get_memsize());
	model->N = 5;
	model->M = 5;
	double vtx[10] = {-4, -1,
			  2.2, -2.6,
			  3.9, -1.8,
			  4, 1,
			  -4.0, 1};
	model->vertex = vtx;
	uint32_t edge[10] = {0, 1,
			     1, 2,
			     2, 3,
			     3, 4,
			     4, 0};
	model->edge = edge;
	
	vcn_mesh_t* mesh = vcn_mesh_create();
	vcn_mesh_generate_from_model(mesh, model);

	vcn_msh3trg_t *msh3trg = vcn_mesh_get_msh3trg(mesh, true, true, true, true, true);

	vcn_mesh_destroy(mesh);

	vcn_msh3trg_save_png(msh3trg, "../../../TRGS.png", 1000, 800);/* TEMPORAL */
	
	CU_ASSERT(false);

	vcn_msh3trg_destroy(msh3trg);
}
