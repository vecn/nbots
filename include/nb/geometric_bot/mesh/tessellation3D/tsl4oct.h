#ifndef __NB_GEOMETRIC_BOT_MESH_TESSELLATION3D_TSL4OCT_H__
#define __NB_GEOMETRIC_BOT_MESH_TESSELLATION3D_TSL4OCT_H__

#include <stdint.h>

#include "nb/geometric_bot/model/model3D.h"
#include "nb/geometric_bot/mesh/tessellation3D/tsl4oct.h"

uint32_t nb_tsl4oct_get_memsize(void);
void nb_tsl4oct_init(void *t3D);
void nb_tsl4oct_copy(void *t3D, const void *src);
void nb_tsl4oct_finish(void *t3D);
void nb_tsl4oct_clear(void* t3D);
void nb_tsl4oct_set_density(void* t3D,
			    double (*density)(const double x[3],
					      const void *data),
			    const void *density_data);
void nb_tsl4oct_unset_density(void* t3D);
bool nb_tsl4oct_is_empty(const void *const t3D);
void nb_tsl4oct_generate_from_model(void *t3D,
				    const nb_model3D_t *const model);
void nb_tsl4oct_get_simplest_from_model(void *t3D,
					const nb_model3D_t *const  model);
bool nb_tsl4oct_is_vtx_inside(const void *const t3D,
			      const double vtx[3]);
void nb_tsl4oct_refine(void *t3D);
bool nb_tsl4oct_insert_vtx(void *t3D,
			   const double vertex[3]);
uint32_t nb_tsl4oct_get_N_vtx(const void *t3D);
uint32_t nb_tsl4oct_get_N_edges(const void *t3D);
uint32_t nb_tsl4oct_get_N_faces(const void *t3D);
uint32_t nb_tsl4oct_get_N_tetra(const void *t3D);
double nb_tsl4oct_get_vol(const void *t3D);

#endif
