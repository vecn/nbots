#include <stdbool.h>
#include <stdlib.h>
#include <alloca.h>

#include <CUnit/Basic.h>

#include "nb/geometric_bot.h"

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
	vcn_model_t *model = vcn_model_create_polygon(1, 0, 0, 6);
	vcn_mesh_t* mesh = vcn_mesh_create();
	vcn_mesh_set_geometric_constraint(mesh,
					 NB_MESH_GEOM_CONSTRAINT_MAX_EDGE_LENGTH,
					 1.5);
	vcn_mesh_generate_from_model(mesh, model);
	vcn_model_destroy(model);
	
	uint32_t size = nb_mshquad_get_memsize();
	nb_mshquad_t *quad = alloca(size);
	nb_mshquad_init(quad);
	nb_mshquad_load_from_mesh(quad, mesh);
	vcn_mesh_destroy(mesh);
	
	CU_ASSERT(7 == quad->N_nod);
	CU_ASSERT(9 == quad->N_edg);
	CU_ASSERT(3 == quad->N_elems);

	nb_mshquad_finish(quad);
}
