#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>
#include <alloca.h>

#include "nb/geometric_bot/mesh/mesh2D.h"
#include "nb/graph_bot.h"

#include "nb/geometric_bot/mesh/elements2D/quad.h"

#include "../mesh2D_structs.h"

typedef struct {
	uint32_t N_matchs;
	uint32_t N_unmatched_trg;
} match_data;

static void set_arrays_memory(nb_mshquad_t *quad);
static uint32_t get_size_of_nod_x_sgm(const nb_mshquad_t *const quad);
static void set_size_of_nod_x_sgm(nb_mshquad_t *quad, char *memblock);
static void copy_nodes(nb_mshquad_t* quad, const nb_mshquad_t *const src_quad);
static void copy_edges(nb_mshquad_t* quad, const nb_mshquad_t *const src_quad);
static void copy_elems(nb_mshquad_t* quad, const nb_mshquad_t *const src_quad);
static void copy_vtx(nb_mshquad_t* quad, const nb_mshquad_t *const src_quad);
static void copy_N_nod_x_sgm(nb_mshquad_t* quad, const nb_mshquad_t *const src_quad);
static void copy_nod_x_sgm(nb_mshquad_t* quad, const nb_mshquad_t *const src_quad);
static void* malloc_quad(void);
static void set_quad_quality_as_weights(const nb_mesh_t *const mesh,
					nb_graph_t *graph);

static void set_mshquad(nb_mshquad_t *quad,
			const nb_graph_t *const graph,
			const nb_mesh_t *const mesh,
			const uint32_t *const matches);
static void get_match_data(const nb_graph_t *const graph,
			   const uint32_t *const matches,
			   match_data *data);
static void set_nodes(nb_mshquad_t *quad, const nb_mesh_t *const mesh);
static void set_edges(nb_mshquad_t *quad, const nb_mesh_t *const mesh,
		      const uint32_t *const matches);
static bool edge_is_not_matched(const msh_edge_t *const edge,
				const uint32_t *const matches);
static void set_elems(nb_mshquad_t *quad, const nb_mesh_t *const mesh,
		      const uint32_t *const matches);
static void set_vtx(nb_mshquad_t *quad, const nb_mesh_t *const mesh);
static void set_N_nod_x_sgm(nb_mshquad_t *quad, const nb_mesh_t *const mesh);
static void set_nod_x_sgm(nb_mshquad_t *quad, const nb_mesh_t *const mesh);

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
	
	set_arrays_memory(quad);

	copy_nodes(quad, src_quad);
	copy_edges(quad, src_quad);
	copy_elems(quad, src_quad);
	copy_vtx(quad, src_quad);
	copy_N_nod_x_sgm(quad, src_quad);

	uint32_t nod_x_sgm_size = get_size_of_nod_x_sgm(quad);
	char *memblock = malloc(nod_x_sgm_size);
	set_size_of_nod_x_sgm(quad, memblock);

	copy_nod_x_sgm(quad, src_quad);
}

static void set_arrays_memory(nb_mshquad_t *quad)
{
	uint32_t nod_size = quad->N_nod * 2 * sizeof(*(quad->nod));
	uint32_t edg_size = quad->N_edg * 2 * sizeof(*(quad->edg));
	uint32_t type_size = quad->N_elems * sizeof(*(quad->type));
	uint32_t adj_size = quad->N_elems * 4 * sizeof(*(quad->adj));
	uint32_t ngb_size = quad->N_elems * 4 * sizeof(*(quad->ngb));
	uint32_t vtx_size = quad->N_vtx * sizeof(*(quad->vtx));
	uint32_t sgm_size = quad->N_sgm * sizeof(*(quad->N_nod_x_sgm));

	uint32_t size = nod_size + edg_size + type_size +
		adj_size + ngb_size + vtx_size + sgm_size;

	char *memblock = malloc(size);

	quad->nod = (void*) memblock;
	quad->edg = (void*) ((char*)(quad->nod) + nod_size);
	quad->type = (void*) ((char*)(quad->edg) + edg_size);
	quad->adj = (void*) ((char*)(quad->type) + type_size);
	quad->ngb = (void*) ((char*)(quad->adj) + adj_size);
	quad->vtx = (void*) ((char*)(quad->ngb) + ngb_size);
	quad->N_nod_x_sgm = (void*) ((char*)(quad->vtx) + vtx_size);
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

static void copy_N_nod_x_sgm(nb_mshquad_t* quad, const nb_mshquad_t *const src_quad)
{
	memcpy(quad->N_nod_x_sgm, src_quad->N_nod_x_sgm,
	       quad->N_sgm * sizeof(*(quad->N_nod_x_sgm)));
}

static void copy_nod_x_sgm(nb_mshquad_t* quad, const nb_mshquad_t *const src_quad)
{	
	for (int i = 0; i < quad->N_sgm; i++) {
		memcpy(&(quad->nod_x_sgm[i]), &(src_quad->nod_x_sgm[i]),
		       quad->N_nod_x_sgm[i] * sizeof(*(quad->nod_x_sgm[i])));
	}

}

void nb_mshquad_finish(void *mshquad_ptr)
{
	nb_mshquad_t *quad = mshquad_ptr;
	if (NULL != quad->nod) {
		free(quad->nod_x_sgm);
		free(quad->nod);		
	}
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
	nb_mshquad_finish(mshquad_ptr);
	free(mshquad_ptr);
}

void nb_mshquad_clear(void *mshquad_ptr)
{
	nb_mshquad_t *quad = mshquad_ptr;
	if (NULL != quad->nod) {
		free(quad->nod_x_sgm);
		free(quad->nod);		
	}
	memset(0, mshquad_ptr, nb_mshquad_get_memsize());	
}

void nb_mshquad_load_from_mesh(nb_mshquad_t *mshquad,
			       const nb_mesh_t *const mesh)
{
	mesh_alloc_vtx_ids((vcn_mesh_t*)mesh);
	mesh_alloc_trg_ids((vcn_mesh_t*)mesh);
	nb_graph_t *graph = vcn_mesh_create_elem_graph(mesh);

	set_quad_quality_as_weights(mesh, graph);

	uint32_t *matches = malloc(graph->N);
	nb_graph_matching_greedy(graph, matches);
	
	set_mshquad(mshquad, graph, mesh, matches);

	free(matches);
	vcn_graph_destroy(graph);
	mesh_free_vtx_ids((vcn_mesh_t*)mesh);
	mesh_free_trg_ids((vcn_mesh_t*)mesh);
}

static void set_quad_quality_as_weights(const nb_mesh_t *const mesh,
					nb_graph_t *graph)
{
	
}

static void set_mshquad(nb_mshquad_t *quad,
			const nb_graph_t *const graph,
			const nb_mesh_t *const mesh,
			const uint32_t *const matches)
{
	match_data *data = alloca(sizeof(match_data));
	get_match_data(graph, matches, data);
	
	quad->N_nod = vcn_mesh_get_N_vtx(mesh);
	quad->N_edg = nb_graph_get_N_edges(graph) - data->N_matchs;
	quad->N_elems = data->N_matchs + data->N_unmatched_trg;
	quad->N_vtx = mesh->N_input_vtx;
	quad->N_sgm = mesh->N_input_sgm;
	
	set_arrays_memory(quad);

	set_nodes(quad, mesh);
	set_edges(quad, mesh, matches);
	set_elems(quad, mesh, matches);
	set_vtx(quad, mesh);

	uint32_t nod_x_sgm_memsize = set_N_nod_x_sgm(quad, mesh);

	char *memblock = malloc(nod_x_sgm_memsize);
	set_size_of_nod_x_sgm(quad, memblock);
	set_nod_x_sgm(quad, mesh);
}

static void set_nodes(nb_mshquad_t *quad, const nb_mesh_t *const mesh)
{
	vcn_bins2D_iter_t* iter = vcn_bins2D_iter_create();
	vcn_bins2D_iter_set_bins(iter, mesh->ug_vtx);
	while (vcn_bins2D_iter_has_more(iter)) {
		const msh_vtx_t* vtx = vcn_bins2D_iter_get_next(iter);
		uint32_t id = *(uint32_t*)((void**)vtx->attr)[0];
		quad->nod[id * 2] = vtx->x[0] / mesh->scale + mesh->xdisp;
		quad->nod[id*2+1] = vtx->x[1] / mesh->scale + mesh->ydisp;
	}
	vcn_bins2D_iter_destroy(iter);	
}

static void set_edges(nb_mshquad_t *quad, const nb_mesh_t *const mesh,
		      const uint32_t *const matches)
{
	uint32_t i = 0;
	uint16_t iter_size = nb_iterator_get_memsize();
	nb_iterator_t *iter = alloca(iter_size);
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, mesh->ht_edge);
	while (nb_iterator_has_more(iter)) {
		msh_edge_t *edge = (msh_edge_t*) nb_iterator_get_next(iter);
		if (edge_is_not_matched(edge, matches)) {
			quad->edg[i * 2] = *(uint32_t*)((void**)edge->v1->attr)[0];
			quad->edg[i*2+1] = *(uint32_t*)((void**)edge->v2->attr)[0];
			i += 1;
		}
	}
	nb_iterator_finish(iter);
}

static bool edge_is_not_matched(const msh_edge_t *const edge,
				const uint32_t *const matches)
{
	int id1 = *(uint32_t*)((void**)edge->t1->attr)[0];
	int id2 = *(uint32_t*)((void**)edge->t2->attr)[0];
	bool out;
	if (matches[id1] != id2)
		out = true;
	else
		out = false;
	return out;
}

static void set_elems(nb_mshquad_t *quad, const nb_mesh_t *const mesh,
		      const uint32_t *const matches)
{
/* AQUI VOY */
	uint16_t iter_size = nb_iterator_get_memsize();
	nb_iterator_t *iter = alloca(iter_size);
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, mesh->ht_trg);
	while (nb_iterator_has_more(iter)) {
		msh_trg_t* trg = nb_iterator_get_next(iter);
		uint32_t id = ((uint32_t*)((void**)trg->attr)[0])[0];

		if (matches[]){}

		uint32_t elem_id = msh3trg->N_triangles;
		if (NULL != trg->t1)
			id1 = ((uint32_t*)((void**)trg->t1->attr)[0])[0];

		uint32_t id2 = msh3trg->N_triangles;
		if (NULL != trg->t2)
			id2 = ((uint32_t*)((void**)trg->t2->attr)[0])[0];
		uint32_t id3 = msh3trg->N_triangles;
		if (NULL != trg->t3)
			id3 = ((uint32_t*)((void**)trg->t3->attr)[0])[0];

		quad->ngb[id * 3] = id1;
		quad->ngb[id*3+1] = id2;
		quad->ngb[id*3+2] = id3;
	}
	nb_iterator_finish(iter);
}

static void set_vtx(nb_mshquad_t *quad, const nb_mesh_t *const mesh);
static void set_N_nod_x_sgm(nb_mshquad_t *quad, const nb_mesh_t *const mesh);
static void set_nod_x_sgm(nb_mshquad_t *quad, const nb_mesh_t *const mesh);

static void get_match_data(const nb_graph_t *const graph,
			   const uint32_t *const matches,
			   match_data *data)
{
	data->N_matchs = 0;
	for (uint32_t i = 0; i < graph->N; i++) {
		if (matches[i] != i)
			data->N_matchs += 1;
		else
			data->N_unmatched_trg += 1;
	}
	data->N_matchs /= 2; /* Double counting */
}

void nb_mshquad_set_nodal_graph(const nb_mshquad_t *mshquad,
				nb_graph_t *graph);
void nb_mshquad_set_elemental_graph(const nb_mshquad_t *mshquad,
				    nb_graph_t *graph);
