#ifndef __NB_GEOMETRIC_BOT_MESH_ELEMENTS2D_QUAD_H__
#define __NB_GEOMETRIC_BOT_MESH_ELEMENTS2D_QUAD_H__

#include <stdint.h>
#include "quad_struct.h"
#include "nb/geometric_bot/mesh/mesh2D.h"
#include "nb/graph_bot.h"

uint32_t nb_quad_get_memsize(void);
void nb_mshquad_init(void *mshquad_ptr);
void nb_mshquad_copy(void *dest, const void *const src);
void nb_mshquad_finish(void *mshquad_ptr);

void* nb_mshquad_create(void);
void* nb_mshquad_clone(const void *const mshquad_ptr);
void nb_mshquad_destroy(void *mshquad_ptr);  

void nb_mshquad_clear(void *mshquad_ptr);

void nb_mshquad_load_from_mesh(nb_mshquad_t *mshquad,
			       const nb_mesh_t *const mesh);

void nb_mshquad_set_nodal_graph(const nb_mshquad_t *mshquad,
				nb_graph_t *graph);
void nb_mshquad_set_elemental_graph(const nb_mshquad_t *mshquad,
				    nb_graph_t *graph);

#endif
