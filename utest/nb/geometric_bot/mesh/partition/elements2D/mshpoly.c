#include <stdbool.h>
#include <stdlib.h>

#include <CUnit/Basic.h>

#include "nb/memory_bot.h"
#include "nb/geometric_bot.h"

static int suite_init(void);
static int suite_clean(void);

static void test_load_from_mesh(void);

void cunit_nb_mshpoly(void)
{
	CU_pSuite suite =
		CU_add_suite("nb/geometric_bot/mesh/partition/elements2D/mshpoly.c",
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
	vcn_model_t *model = vcn_model_create_polygon(1, 0, 0, 6);
	nb_mesh_t* mesh = nb_mesh_create();
	nb_mesh_set_geometric_constraint(mesh,
					  NB_MESH_GEOM_CONSTRAINT_MAX_EDGE_LENGTH,
					  0.05);
	nb_mesh_generate_from_model(mesh, model);
	vcn_model_destroy(model);

	void *poly = nb_allocate_on_stack(nb_mshpoly_get_memsize());
	nb_mshpoly_init(poly);
	nb_mshpoly_load_from_mesh(poly, mesh);
	nb_mesh_destroy(mesh);

	/* TEMPORAL FAIL: Produce different triangles each time */
	uint32_t N_elems = nb_mshpoly_get_N_elems(poly);
	uint32_t N_edges = nb_mshpoly_get_N_edges(poly);
	CU_ASSERT(2000 < N_elems && 2900 > N_elems);
	CU_ASSERT(7500 < N_edges && 8300 > N_edges);

	nb_mshpoly_finish(poly);
}
