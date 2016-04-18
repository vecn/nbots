/*
 * Compile with:
gcc -g -std=c99 -I ../include -L libs/nbots/shared/debug/ -L libs/nbots_cairo/shared/debug nbots_tmp.c -lm -lnbots -lnbots_cairo
export LD_LIBRARY_PATH=libs/nbots/shared/debug:libs/nbots_cairo/shared/debug
./a.out
 */

#include <stdbool.h>
#include <stdlib.h>
#include <alloca.h>


#include "nb/geometric_bot/point2D.h"

#define OUTPUT_DIR "./"

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

static void test_load_from_mesh(void);

int main(void)
{
	test_load_from_mesh();
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
}
