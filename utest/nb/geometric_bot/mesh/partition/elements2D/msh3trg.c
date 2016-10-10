#include <stdbool.h>
#include <stdlib.h>

#include <CUnit/Basic.h>

#include "nb/memory_bot.h"
#include "nb/graphics_bot.h"
#include "nb/geometric_bot.h"

static int suite_init(void);
static int suite_clean(void);

static void test_load_from_mesh(void);

void cunit_nb_msh3trg(void)
{
	CU_pSuite suite = CU_add_suite("nb/geometric_bot/mesh/partition/elements2D/msh3trg.c",
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
	nb_model_t *model  = nb_allocate_on_stack(nb_model_get_memsize());
	model->N = 5;
	model->M = 5;
	model->H = 0;
	double vtx[10] = {-4, -1,
			  2.2, -2.6,
			  3.9, -2.7,
			  4, 1,
			  -4.0, 1};
	model->vertex = vtx;
	uint32_t edge[10] = {0, 1,
			     1, 2,
			     2, 3,
			     3, 4,
			     4, 0};
	model->edge = edge;
	model->holes = NULL;
	
	uint32_t mesh_memsize = nb_mesh_get_memsize();
	nb_mesh_t* mesh = nb_allocate_on_stack(mesh_memsize);
	nb_mesh_init(mesh);
	nb_mesh_generate_from_model(mesh, model);

	uint32_t msh_memsize = nb_msh3trg_get_memsize();
	void *msh3trg = nb_allocate_on_stack(msh_memsize);
	nb_msh3trg_init(msh3trg);
	nb_msh3trg_load_from_mesh(msh3trg, mesh);

	nb_mesh_finish(mesh);
	
	CU_ASSERT(9 == nb_msh3trg_get_N_nodes(msh3trg));
	CU_ASSERT(15 == nb_msh3trg_get_N_edges(msh3trg));
	CU_ASSERT(7 == nb_msh3trg_get_N_elems(msh3trg));
	
	nb_msh3trg_finish(msh3trg);
}
