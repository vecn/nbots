#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "nb/geometric_bot/mesh/mesh2D.h"
#include "nb/graph_bot.h"

#include "nb/geometric_bot/mesh/elements2D/quad.h"

static void set_arrays_memory(nb_mshquad_t *quad,
			      const nb_mshquad_t *const src_quad);
static uint32_t get_size_of_nod_x_sgm(const nb_mshquad_t *const quad);
static void set_size_of_nod_x_sgm(nb_mshquad_t *quad, char *memblock);
static void copy_nodes(nb_mshquad_t* quad, const nb_mshquad_t *const src_quad);
static void copy_edges(nb_mshquad_t* quad, const nb_mshquad_t *const src_quad);
static void copy_elems(nb_mshquad_t* quad, const nb_mshquad_t *const src_quad);
static void copy_vtx(nb_mshquad_t* quad, const nb_mshquad_t *const src_quad);
static void copy_sgm(nb_mshquad_t* quad, const nb_mshquad_t *const src_quad);
static void* malloc_quad(void);

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
	memcpy(dest, src, nb_mshquad_get_memsize());
	nb_mshquad_t *quad = dest;
	const nb_mshquad_t *const src_quad = src;
	
	set_arrays_memory(quad, src_quad);

	copy_nodes(quad, src_quad);
	copy_edges(quad, src_quad);
	copy_elems(quad, src_quad);
	copy_vtx(quad, src_quad);
	copy_sgm(quad, src_quad);
}

static void set_arrays_memory(nb_mshquad_t *quad,
			      const nb_mshquad_t *const src_quad)
{
	uint32_t nod_size = quad->N_nod * 2 * sizeof(*(quad->nod));
	uint32_t edg_size = quad->N_edg * 2 * sizeof(*(quad->edg));
	uint32_t type_size = quad->N_elems * sizeof(*(quad->type));
	uint32_t adj_size = quad->N_elems * 4 * sizeof(*(quad->adj));
	uint32_t ngb_size = quad->N_elems * 4 * sizeof(*(quad->ngb));
	uint32_t vtx_size = quad->N_vtx * sizeof(*(quad->vtx));
	uint32_t sgm_size = quad->N_sgm * sizeof(*(quad->N_nod_x_sgm));

	uint32_t size = nod_size + edg_size + type_size +
		adj_size + ngb_size + vtx_size +
		sgm_size + get_size_of_nod_x_sgm(quad);

	char *memblock = malloc(size);

	quad->nod = (void*) memblock;
	quad->edg = (void*) ((char*)(quad->nod) + nod_size);
	quad->type = (void*) ((char*)(quad->edg) + edg_size);
	quad->adj = (void*) ((char*)(quad->type) + type_size);
	quad->ngb = (void*) ((char*)(quad->adj) + adj_size);
	quad->vtx = (void*) ((char*)(quad->ngb) + ngb_size);
	quad->N_nod_x_sgm = (void*) ((char*)(quad->vtx) + vtx_size);
	set_size_of_nod_x_sgm(quad, (char*)(quad->N_nod_x_sgm) + sgm_size);
}

static uint32_t get_size_of_nod_x_sgm(const nb_mshquad_t *const quad)
{
	uint32_t size = 0;
	for (uint32_t i = 0; i < quad->N_sgm; i++)
		size += quad->N_nod_x_sgm[i] *
			sizeof(*(quad->nod_x_sgm[i]));
	return size;
}

static void set_size_of_nod_x_sgm(nb_mshquad_t *quad, char *memblock)
{
	for (uint32_t i = 0; i < quad->N_sgm; i++) {
		quad->nod_x_sgm[i] = memblock;
		memblock += quad->N_nod_x_sgm[i] *
			sizeof(*(quad->nod_x_sgm[i]));
	}
}

static void copy_nodes(nb_mshquad_t* quad, const nb_mshquad_t *const src_quad)
{
	memcpy(quad->nod, src_quad->nod,
	       2 * quad->N_nod * sizeof(*(quad->nod)));
}

static void copy_edges(nb_mshquad_t* quad, const nb_mshquad_t *const src_quad)
{
	memcpy(quad->edg, src_quad->edg,
	       2 * quad->N_edg * sizeof(*(quad->edg)));
}

static void copy_elems(nb_mshquad_t* quad, const nb_mshquad_t *const src_quad)
{
	memcpy(quad->type, src_quad->type,
	       quad->N_elems * sizeof(*(quad->type)));
	memcpy(quad->adj, src_quad->adj,
	       4 * quad->N_elems * sizeof(*(quad->adj)));
	memcpy(quad->ngb, src_quad->ngb,
	       4 * quad->N_elems * sizeof(*(quad->ngb)));
}

static void copy_vtx(nb_mshquad_t* quad, const nb_mshquad_t *const src_quad)
{
	memcpy(quad->vtx, src_quad->vtx,
	       quad->N_vtx * sizeof(*(quad->vtx)));
}

static void copy_sgm(nb_mshquad_t* quad, const nb_mshquad_t *const src_quad)
{
	memcpy(quad->N_nod_x_sgm, src_quad->N_nod_x_sgm,
	       quad->N_sgm * sizeof(*(quad->N_nod_x_sgm)));
	
	for (int i = 0; i < quad->N_sgm; i++) {
		memcpy(&(quad->nod_x_sgm[i]), &(src_quad->nod_x_sgm[i]),
		       quad->N_nod_x_sgm[i] * sizeof(*(quad->nod_x_sgm[i])));
	}

}

void nb_mshquad_clear(void *mshquad_ptr)
{
	nb_mshquad_t *quad = mshquad_ptr;
	if (NULL != quad->nod)
		free(quad->nod);
	memset(0, mshquad_ptr, nb_mshquad_get_memsize());	
}

void* nb_mshquad_create(void)
{
	nb_mshquad_t *quad = malloc_mshquad();
	nb_mshquad_init(quad);
	return quad;
}

static void* malloc_quad(void)
{
	uint32_t size = nb_mshquad_get_memsize();
	nb_mshquad_t *quad = malloc(size);
	return quad;
}

void* nb_mshquad_clone(const void *const mshquad_ptr)
{
	nb_mshquad_t *quad = malloc_mshquad();
	nb_mshquad_copy(quad, mshquad_ptr);
	return quad;
}

void nb_mshquad_destroy(void *mshquad_ptr)
{
	nb_mshquad_clear(mshquad_ptr);
	free(mshquad_ptr);
}

void nb_mshquad_load_from_mesh(nb_mshquad_t *mshquad,
			       const nb_mesh_t *const mesh)
{
	nb_graph_t *graph = vcn_mesh_create_elem_graph(mesh);

	set_quad_quality_as_weights(mesh, graph);

	int8_t *edges = malloc(graph->N_adj);
	nb_graph_matching_greedy(graph, edges);
	
	/* AQUI VOY assemble_mshquad(graph, mesh, edges)*/

	free(edges);
}

void nb_mshquad_set_nodal_graph(const nb_mshquad_t *mshquad,
				nb_graph_t *graph);
void nb_mshquad_set_elemental_graph(const nb_mshquad_t *mshquad,
				    nb_graph_t *graph);
