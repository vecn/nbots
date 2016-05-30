#ifndef __NB_GEOMETRIC_BOT_MESH_ELEMENTS2D_POLYGONS_H__
#define __NB_GEOMETRIC_BOT_MESH_ELEMENTS2D_POLYGONS_H__

#include <stdbool.h>
#include <stdint.h>

#include "polygons_struct.h"
#include "nb/geometric_bot/mesh/mesh2D.h"
#include "nb/graph_bot.h"

uint32_t nb_mshpoly_get_memsize(void);
void nb_mshpoly_init(void *mshpoly_ptr);
void nb_mshpoly_copy(void *dest, const void *const src);
void nb_mshpoly_finish(void *mshpoly_ptr);

void* nb_mshpoly_create(void);
void* nb_mshpoly_clone(const void *const mshpoly_ptr);
void nb_mshpoly_destroy(void *mshpoly_ptr);  

void nb_mshpoly_clear(void *mshpoly_ptr);

void nb_mshpoly_load_from_mesh(nb_mshpoly_t *mshpoly, nb_mesh_t *mesh);

void nb_mshpoly_Lloyd_iteration(nb_mshpoly_t *mshpoly, uint32_t max_iter,
				double (*density)(const double[2],
						  const void *data),
				const void *density_data);

void nb_mshpoly_set_fem_graph(const nb_mshpoly_t *mshpoly,
				nb_graph_t *graph);
void nb_mshpoly_set_nodal_graph(const nb_mshpoly_t *mshpoly,
				nb_graph_t *graph);
void nb_mshpoly_set_elemental_graph(const nb_mshpoly_t *mshpoly,
				    nb_graph_t *graph);


#endif
