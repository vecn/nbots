#include <stdint.h>
#include <stdbool.h>

#include "nb/geometric_bot/model/model3D.h"
#include "nb/geometric_bot/mesh/tessellator3D.h"
#include "nb/geometric_bot/mesh/tessellation3D/tsl4oct.h"

#include "tessellator3D_struct.h"

static void set_tessellator_interface(nb_tessellator3D_t *tessel3D,
				      nb_tessellator3D_type  type);
static void set_tsl4oct_interface(nb_tessellator3D_t *tessel3D);

uint32_t nb_tessellator3D_get_memsize(nb_tessellator3D_type type)
{
	uint32_t mem;
	switch (type) {
	case NB_OCTREE:
		mem = nb_tsl4oct_get_memsize();
		break;
	default:
		mem = nb_tsl4oct_get_memsize();
		break;
	}
	return mem + sizeof(nb_tessellator3D_t);
}

void nb_tessellator3D_init(nb_tessellator3D_t *tessel3D,
			   nb_tessellator3D_type type)
{
	char *memblock = (void*) tessel3D;
	tessel3D->t3D = (void*) (memblock + sizeof(nb_tessellator3D_t));
	tessel3D->type = type;
	set_tessellator_interface(tessel3D, type);
	tessel3D->init(tessel3D->t3D);
}

static void set_tessellator_interface(nb_tessellator3D_t *tessel3D,
				      nb_tessellator3D_type  type)
{
	switch (type) {
	case NB_OCTREE:
		set_tsl4oct_interface(tessel3D);
		break;
	default:
		set_tsl4oct_interface(tessel3D);
		break;
	}
}

static void set_tsl4oct_interface(nb_tessellator3D_t *tessel3D)
{
	tessel3D->get_memsize = nb_tsl4oct_get_memsize;
	tessel3D->init = nb_tsl4oct_init;
	tessel3D->copy = nb_tsl4oct_copy;
	tessel3D->finish = nb_tsl4oct_finish;
	tessel3D->clear = nb_tsl4oct_clear;
	tessel3D->set_density = nb_tsl4oct_set_density;
	tessel3D->unset_density = nb_tsl4oct_unset_density;
	tessel3D->is_empty = nb_tsl4oct_is_empty;
	tessel3D->generate_from_model = nb_tsl4oct_generate_from_model;
	tessel3D->get_simplest_from_model = nb_tsl4oct_get_simplest_from_model;
	tessel3D->is_vtx_inside = nb_tsl4oct_is_vtx_inside;
	tessel3D->refine = nb_tsl4oct_refine;
	tessel3D->insert_vtx = nb_tsl4oct_insert_vtx;
	tessel3D->get_N_vtx = nb_tsl4oct_get_N_vtx;
	tessel3D->get_N_edges = nb_tsl4oct_get_N_edges;
	tessel3D->get_N_faces = nb_tsl4oct_get_N_faces;
	tessel3D->get_N_tetra = nb_tsl4oct_get_N_tetra;
	tessel3D->get_vol = nb_tsl4oct_get_vol;
}

void nb_tessellator3D_copy(nb_tessellator3D_t *tessel3D, 
			  const nb_tessellator3D_t *src)
{
	tessel3D->copy(tessel3D->t3D, src->t3D);
}

void nb_tessellator3D_finish(nb_tessellator3D_t *tessel3D)
{
	tessel3D->finish(tessel3D->t3D);
}

void nb_tessellator3D_clear(nb_tessellator3D_t* tessel3D)
{
	tessel3D->clear(tessel3D->t3D);
}

void nb_tessellator3D_set_density(nb_tessellator3D_t* tessel3D,
				  double (*density)(const double x[3],
						    const void *data),
				  const void *density_data)
{
	tessel3D->set_density(tessel3D->t3D, density, density_data);
}

void nb_tessellator3D_unset_density(nb_tessellator3D_t* tessel3D)
{
	tessel3D->unset_density(tessel3D->t3D);
}

bool nb_tessellator3D_is_empty(const nb_tessellator3D_t *const tessel3D)
{
	return tessel3D->is_empty(tessel3D->t3D);
}

void nb_tessellator3D_generate_from_model(nb_tessellator3D_t *tessel3D,
					  const nb_model3D_t *const model)
{
	tessel3D->generate_from_model(tessel3D->t3D, model);
}

void nb_tessellator3D_get_simplest_from_model(nb_tessellator3D_t *tessel3D,
					      const nb_model3D_t *const  model)
{
	tessel3D->get_simplest_from_model(tessel3D->t3D, model);
}

bool nb_tessellator3D_is_vtx_inside(const nb_tessellator3D_t *const tessel3D,
				    const double vtx[3])
{
	return tessel3D->is_vtx_inside(tessel3D->t3D, vtx);
}

void nb_tessellator3D_refine(nb_tessellator3D_t *tessel3D)
{
	tessel3D->refine(tessel3D->t3D);
}

bool nb_tessellator3D_insert_vtx(nb_tessellator3D_t *tessel3D,
				 const double vertex[3])
{
	return tessel3D->insert_vtx(tessel3D->t3D, vertex);
}

uint32_t nb_tessellator3D_get_N_vtx(const nb_tessellator3D_t *tessel3D)
{
	return tessel3D->get_N_vtx(tessel3D->t3D);
}

uint32_t nb_tessellator3D_get_N_edges(const nb_tessellator3D_t *tessel3D)
{
	return tessel3D->get_N_edges(tessel3D->t3D);
}

uint32_t nb_tessellator3D_get_N_faces(const nb_tessellator3D_t *tessel3D)
{
	return tessel3D->get_N_faces(tessel3D->t3D);
}

uint32_t nb_tessellator3D_get_N_tetra(const nb_tessellator3D_t *tessel3D)
{
	return tessel3D->get_N_tetra(tessel3D->t3D);
}

double nb_tessellator3D_get_vol(const nb_tessellator3D_t *tessel3D)
{
	return tessel3D->get_vol(tessel3D->t3D);
}
