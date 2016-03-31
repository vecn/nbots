#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "nb/geometric_bot/mesh/mesh2D.h"
#include "nb/graph_bot.h"

#include "nb/geometric_bot/mesh/elements2D/quad.h"

static void copy_nodes(nb_mshquad_t* quad, const nb_mshquad_t *const src_quad);
static void copy_edges(nb_mshquad_t* quad, const nb_mshquad_t *const src_quad);
static void copy_elems(nb_mshquad_t* quad, const nb_mshquad_t *const src_quad);
static void copy_vtx(nb_mshquad_t* quad, const nb_mshquad_t *const src_quad);
static void copy_sgm(nb_mshquad_t* quad, const nb_mshquad_t *const src_quad);

uint32_t nb_mshquad_get_memsize(void)
{
	return sizeof(nb_mshquad_t);
}

void nb_mshquad_init(void *mshquad_ptr)
{
	memset(0, mshquad_ptr, nb_mshquad_get_memsize());
}

void nb_mshquad_copy(void *dest, const void *const src)
{
	nb_mshquad_t *quad = dest;
	const nb_mshquad_t *const src_quad = src;
	memcpy(dest, src, nb_mshquad_get_memsize());
	
	copy_nodes(quad, src_quad);
	copy_edges(quad, src_quad);
	copy_elems(quad, src_quad);
	copy_vtx(quad, src_quad);
	copy_sgm(quad, src_quad);
}

static void copy_nodes(nb_mshquad_t* quad, const nb_mshquad_t *const src_quad)
{
	memcpy(quad->nod, src_quad->nod,
	       2 * quad->N_nod * sizeof(*(quad->nod)));
}

static void copy_edges(nb_mshquad_t* quad, const nb_mshquad_t *const src_quad);
static void copy_elems(nb_mshquad_t* quad, const nb_mshquad_t *const src_quad);
static void copy_vtx(nb_mshquad_t* quad, const nb_mshquad_t *const src_quad);
static void copy_sgm(nb_mshquad_t* quad, const nb_mshquad_t *const src_quad);

void nb_mshquad_clear(void *mshquad_ptr);

void* nb_mshquad_create(void);
void* nb_mshquad_clone(const void *const mshquad_ptr);
void nb_mshquad(void *mshquad_ptr);  

void nb_mshquad_load_from_mesh(nb_mshquad_t *mshquad,
			       const nb_mesh_t *const mesh);

void nb_mshquad_set_nodal_graph(const nb_mshquad_t *mshquad,
				nb_graph_t *graph);
void nb_mshquad_set_elemental_graph(const nb_mshquad_t *mshquad,
				    nb_graph_t *graph);
