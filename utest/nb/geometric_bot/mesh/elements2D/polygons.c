#include <stdbool.h>
#include <stdlib.h>
#include <alloca.h>

#include <CUnit/Basic.h>

#include "nb/geometric_bot.h"

static int suite_init(void);
static int suite_clean(void);

static void test_load_from_mesh(void);

void cunit_nb_geometric_bot_mesh_elements2D_poly(void)
{
	CU_pSuite suite =
		CU_add_suite("nb/geometric_bot/mesh/elements2D/polygons.c",
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
	vcn_mesh_t* mesh = vcn_mesh_create();
	vcn_mesh_set_geometric_constraint(mesh,
					  NB_MESH_GEOM_CONSTRAINT_MAX_SUBSGM_LENGTH,
					  0.02);
	vcn_mesh_set_geometric_constraint(mesh,
					  NB_MESH_GEOM_CONSTRAINT_MAX_EDGE_LENGTH,
					  0.05);
	vcn_mesh_generate_from_model(mesh, model);
	vcn_model_destroy(model);

	void *poly = alloca(nb_mshpoly_get_memsize());
	nb_mshpoly_init(poly);
	nb_mshpoly_load_from_mesh(poly, mesh);
	vcn_mesh_destroy(mesh);

	/* TEMPORAL FAIL: Produce different triangles each time */
	uint32_t N_elems = nb_mshpoly_get_N_elems(poly);
	uint32_t N_edges = nb_mshpoly_get_N_edges(poly);
	printf("\n --N: %zu %zu\n", N_elems, N_edges);/* TEMPORAL */
	CU_ASSERT(2750 < N_elems && 3150 > N_elems);
	CU_ASSERT(8700 < N_edges && 9200 > N_edges);

	nb_mshpoly_finish(poly);
}
