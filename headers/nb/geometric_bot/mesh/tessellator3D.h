#ifndef __NB_GEOMETRIC_BOT_MESH_TESSELLATOR3D_H__
#define __NB_GEOMETRIC_BOT_MESH_TESSELLATOR3D_H__

#include <stdint.h>
#include <stdbool.h>

#include "nb/geometric_bot/model/model3D.h"

typedef enum {
	NB_DELAUNAY, NB_OCTREE
} nb_tessellator3D_type;

typedef struct nb_tessellator3D_s nb_tessellator3D_t;

uint32_t nb_tessellator3D_get_memsize(nb_tessellator3D_type type);
void nb_tessellator3D_init(nb_tessellator3D_t *tessel3D,
			   nb_tessellator3D_type type);
void nb_tessellator3D_copy(nb_tessellator3D_t *tessel3D, 
			  const nb_tessellator3D_t *src);
void nb_tessellator3D_finish(nb_tessellator3D_t *tessel3D);

void nb_tessellator3D_clear(nb_tessellator3D_t* tessel3D);
void nb_tessellator3D_set_density(nb_tessellator3D_t* tessel3D,
				  double (*density)(const double x[3],
						    const void *data),
				  const void *density_data);
void nb_tessellator3D_unset_density(nb_tessellator3D_t* tessel3D);
bool nb_tessellator3D_is_empty(const nb_tessellator3D_t *const tessel3D);
void nb_tessellator3D_generate_from_model(nb_tessellator3D_t *tessel3D,
					  const nb_model3D_t *const model);

void nb_tessellator3D_get_simplest_from_model(nb_tessellator3D_t *tessel3D,
					      const nb_model3D_t *const  model);
bool nb_tessellator3D_is_vtx_inside(const nb_tessellator3D_t *const tessel3D,
				    const double vtx[3]);
void nb_tessellator3D_refine(nb_tessellator3D_t *tessel3D);
bool nb_tessellator3D_insert_vtx(nb_tessellator3D_t *tessel3D,
				 const double vertex[3]);
uint32_t nb_tessellator3D_get_N_vtx(const nb_tessellator3D_t *tessel3D);
uint32_t nb_tessellator3D_get_N_edges(const nb_tessellator3D_t *tessel3D);
uint32_t nb_tessellator3D_get_N_faces(const nb_tessellator3D_t *tessel3D);
uint32_t nb_tessellator3D_get_N_tetra(const nb_tessellator3D_t *tessel3D);
double nb_tessellator3D_get_vol(const nb_tessellator3D_t *tessel3D);

#endif
