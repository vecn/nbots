#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "nb/geometric_bot/point2D.h"
#include "nb/geometric_bot/utils2D.h"
#include "nb/geometric_bot/knn/bins2D.h"
#include "nb/geometric_bot/knn/bins2D_iterator.h"
#include "nb/geometric_bot/mesh/mesh2D.h"
#include "nb/geometric_bot/mesh/dewall.h"
#include "nb/geometric_bot/mesh/modules2D/graph_generator.h"
#include "nb/geometric_bot/mesh/elements2D/polygons.h"

#include "../mesh2D_structs.h"

#define POW2(a) ((a)*(a))

typedef struct {
	uint32_t N;
	uint32_t *N_adj;
	msh_edge_t **adj;
} vgraph_t;

static void set_arrays_memory(nb_mshpoly_t *poly);
static void copy_nodes(nb_mshpoly_t* poly,
		       const nb_mshpoly_t *const src_poly);
static void copy_edges(nb_mshpoly_t* poly,
		       const nb_mshpoly_t *const src_poly);
static void copy_N_side(nb_mshpoly_t* poly,
			const nb_mshpoly_t *const src_poly);
static uint32_t get_size_of_adj_and_ngb(const nb_mshpoly_t *const poly);
static void set_mem_of_adj_and_ngb(nb_mshpoly_t *poly, uint32_t memsize);
static void copy_adj_and_ngb(nb_mshpoly_t* poly,
			     const nb_mshpoly_t *const src_poly);

static void copy_elem_vtx(nb_mshpoly_t* poly,
			  const nb_mshpoly_t *const src_poly);
static void copy_N_nod_x_sgm(nb_mshpoly_t* poly,
			     const nb_mshpoly_t *const src_poly);
static uint32_t get_size_of_nod_x_sgm(const nb_mshpoly_t *const poly);
static void set_mem_of_nod_x_sgm(nb_mshpoly_t *poly, uint32_t memsize);
static void copy_nod_x_sgm(nb_mshpoly_t* poly,
			   const nb_mshpoly_t *const src_poly);
static void* malloc_poly(void);
static void init_voronoi_graph(vgraph_t *vgraph,
			       const nb_mesh_t *const mesh);
static void count_vgraph_adj(vgraph_t *vgraph, const nb_mesh_t *const mesh);
static void set_vgraph_adj_mem(vgraph_t *vgraph, char *memblock);
static void set_vgraph_adj(vgraph_t *vgraph, const nb_mesh_t *const mesh,
			   char *memblock);
static bool adj_is_cocircular(const msh_edge_t *const edg);
static bool edge_has_trg_with_cl_ccenter(const msh_edge_t *const edg);
static bool edge_is_encroached(const msh_edge_t *edg, const msh_trg_t *trg);
static void finish_voronoi_graph(vgraph_t *vgraph);
static void set_voronoi(nb_mshpoly_t *poly,
			const vgraph_t *const graph);
static void count_voronoi_data(nb_mshpoly_t *poly,
			       const vgraph_t *const vgraph);
static bool cell_is_boundary(const vgraph_t *const graph,
			     uint32_t i);
static void set_nodes(nb_mshpoly_t *poly,
		      const nb_mesh_t *const mesh,
		      const uint32_t *cc_map);
static void set_edges(nb_mshpoly_t *poly,
		      const nb_mesh_t *const mesh);
static uint32_t set_N_sides(nb_mshpoly_t *poly,
			   const nb_mesh_t *const mesh);
static void set_adj_and_ngb(nb_mshpoly_t *poly,
			    const nb_mesh_t *const mesh);
static void set_elem_vtx(nb_mshpoly_t *poly,
			 const nb_mesh_t *const mesh);
static uint32_t set_N_nod_x_sgm(nb_mshpoly_t *poly,
				const nb_mesh_t *const mesh);
static void set_nod_x_sgm(nb_mshpoly_t *poly,
			  const nb_mesh_t *const mesh);

uint32_t nb_mshpoly_get_memsize(void)
{
	return sizeof(nb_mshpoly_t);
}

void nb_mshpoly_init(void *mshpoly_ptr)
{
	memset(mshpoly_ptr, 0, nb_mshpoly_get_memsize());
}

void nb_mshpoly_copy(void *dest, const void *const src)
{
	memcpy(dest, src, nb_mshpoly_get_memsize());
	nb_mshpoly_t *poly = dest;
	const nb_mshpoly_t *const src_poly = src;

	if (dest->N_elems > 0) {
		set_arrays_memory(poly);

		copy_nodes(poly, src_poly);
		copy_edges(poly, src_poly);
		copy_N_side(poly, src_poly);

		uint32_t memsize = get_size_of_adj_and_ngb(poly);
		set_mem_of_adj_and_ngb(poly, memsize);
		copy_adj_and_ngb(poly, src_poly);

		copy_elem_vtx(poly, src_poly);
		copy_N_nod_x_sgm(poly, src_poly);

		memsize = get_size_of_nod_x_sgm(poly);
		set_mem_of_nod_x_sgm(poly, memsize);
		copy_nod_x_sgm(poly, src_poly);
	}
}

static void set_arrays_memory(nb_mshpoly_t *poly)
{
	uint32_t nod_size = poly->N_nod * 2 * sizeof(*(poly->nod));
	uint32_t edg_size = poly->N_edg * 2 * sizeof(*(poly->edg));
	uint32_t N_side_size = poly->N_elems * sizeof(*(poly->N_side));
	uint32_t adj_size = poly->N_elems * sizeof(*(poly->adj));
	uint32_t ngb_size = poly->N_elems * sizeof(*(poly->ngb));
	uint32_t elem_vtx_size = poly->N_vtx * sizeof(*(poly->elem_vtx));
	uint32_t N_nod_x_sgm_size = poly->N_sgm * sizeof(*(poly->N_nod_x_sgm));
	uint32_t nod_x_sgm_size = poly->N_sgm * sizeof(*(poly->nod_x_sgm));

	uint32_t size = nod_size + edg_size + N_side_size +
		adj_size + ngb_size + elem_vtx_size +
		N_nod_x_sgm_size + nod_x_sgm_size;

	char *memblock = malloc(size);

	poly->nod = (void*) memblock;
	poly->edg = (void*) ((char*)(poly->nod) + nod_size);
	poly->N_side = (void*) ((char*)(poly->edg) + edg_size);
	poly->adj = (void*) ((char*)(poly->N_side) + N_side_size);
	poly->ngb = (void*) ((char*)(poly->adj) + adj_size);
	poly->elem_vtx = (void*) ((char*)(poly->ngb) + ngb_size);
	poly->N_nod_x_sgm = (void*) ((char*)(poly->elem_vtx) +
				     elem_vtx_size);
	poly->nod_x_sgm = (void*) ((char*)(poly->N_nod_x_sgm) +
				   N_nod_x_sgm_size);
}

static void copy_nodes(nb_mshpoly_t* poly, const nb_mshpoly_t *const src_poly)
{
	memcpy(poly->nod, src_poly->nod,
	       2 * poly->N_nod * sizeof(*(poly->nod)));
}

static void copy_edges(nb_mshpoly_t* poly, const nb_mshpoly_t *const src_poly)
{
	memcpy(poly->edg, src_poly->edg,
	       2 * poly->N_edg * sizeof(*(poly->edg)));
}

static void copy_N_side(nb_mshpoly_t* poly, const nb_mshpoly_t *const src_poly)
{
	memcpy(poly->N_side, src_poly->N_side,
	       poly->N_elems * sizeof(*(poly->N_side)));
}

static uint32_t get_size_of_adj_and_ngb(const nb_mshpoly_t *const poly)
{
	uint32_t size = 0;
	for (uint32_t i = 0; i < poly->N_sgm; i++)
		size += poly->N_side[i] * sizeof(**(poly->adj));
	return 2 * size;
}

static void set_mem_of_adj_and_ngb(nb_mshpoly_t *poly, uint32_t memsize)
{
	char *memblock1 = malloc(memsize);
	char *memblock2 = memblock1 + memsize / 2;
	for (uint32_t i = 0; i < poly->N_sgm; i++) {
		poly->adj[i] = (void*) memblock1;
		poly->ngb[i] = (void*) memblock2;
		memblock1 += poly->N_side[i] * sizeof(**(poly->adj));
		memblock2 += poly->N_side[i] * sizeof(**(poly->ngb));
	}
}

static void copy_adj_and_ngb(nb_mshpoly_t* poly,
			     const nb_mshpoly_t *const src_poly)
{
	for (int i = 0; i < poly->N_sgm; i++) {
		memcpy(&(poly->adj[i]), &(src_poly->adj[i]),
		       poly->N_side[i] * sizeof(**(poly->adj)));
		memcpy(&(poly->ngb[i]), &(src_poly->ngb[i]),
		       poly->N_side[i] * sizeof(**(poly->ngb)));
	}

}

static void copy_elem_vtx(nb_mshpoly_t* poly,
			  const nb_mshpoly_t *const src_poly)
{
	memcpy(poly->elem_vtx, src_poly->elem_vtx,
	       poly->N_vtx * sizeof(*(poly->elem_vtx)));
}

static void copy_N_nod_x_sgm(nb_mshpoly_t* poly,
			     const nb_mshpoly_t *const src_poly)
{
	memcpy(poly->N_nod_x_sgm, src_poly->N_nod_x_sgm,
	       poly->N_sgm * sizeof(*(poly->N_nod_x_sgm)));
}


static uint32_t get_size_of_nod_x_sgm(const nb_mshpoly_t *const poly)
{
	uint32_t size = 0;
	for (uint32_t i = 0; i < poly->N_sgm; i++)
		size += poly->N_nod_x_sgm[i] *
			sizeof(**(poly->nod_x_sgm));
	return size;
}

static void set_mem_of_nod_x_sgm(nb_mshpoly_t *poly, uint32_t memsize)
{
	char *memblock = malloc(memsize);
	for (uint32_t i = 0; i < poly->N_sgm; i++) {
		poly->nod_x_sgm[i] = (void*) memblock;
		memblock += poly->N_nod_x_sgm[i] *
			sizeof(**(poly->nod_x_sgm));
	}
}

static void copy_nod_x_sgm(nb_mshpoly_t* poly,
			   const nb_mshpoly_t *const src_poly)
{
	for (int i = 0; i < poly->N_sgm; i++) {
		memcpy(&(poly->nod_x_sgm[i]), &(src_poly->nod_x_sgm[i]),
		       poly->N_nod_x_sgm[i] * sizeof(**(poly->nod_x_sgm)));
	}
}

void nb_mshpoly_finish(void *mshpoly_ptr)
{
	nb_mshpoly_clear(mshpoly_ptr);
}

void* nb_mshpoly_create(void)
{
	nb_mshpoly_t *poly = malloc_poly();
	nb_mshpoly_init(poly);
	return poly;
}

static void* malloc_poly(void)
{
	uint32_t size = nb_mshpoly_get_memsize();
	nb_mshpoly_t *poly = malloc(size);
	return poly;
}

void* nb_mshpoly_clone(const void *const mshpoly_ptr)
{
	nb_mshpoly_t *poly = malloc_poly();
	nb_mshpoly_copy(poly, mshpoly_ptr);
	return poly;
}

void nb_mshpoly_destroy(void *mshpoly_ptr)
{
	nb_mshpoly_finish(mshpoly_ptr);
	free(mshpoly_ptr);
}

void nb_mshpoly_clear(void *mshpoly_ptr)
{
	nb_mshpoly_t *poly = mshpoly_ptr;
	if (NULL != poly->nod) {
		free(poly->adj[0]);
		free(poly->nod_x_sgm[0]);
		free(poly->nod);		
	}
	memset(mshpoly_ptr, 0, nb_mshpoly_get_memsize());
}

void nb_mshpoly_load_from_mesh(nb_mshpoly_t *mshpoly,
			       const nb_mesh_t *const mesh)
{
	if (vcn_mesh_get_N_trg(mesh) > 0) {
		mesh_alloc_vtx_ids((vcn_mesh_t*)mesh);
	
		vgraph_t vgraph;
		init_voronoi_graph(&vgraph, mesh);

		set_voronoi(mshpoly, &vgraph);

		finish_voronoi_graph(&vgraph);

		mesh_free_vtx_ids((vcn_mesh_t*)mesh);
	}
}

static void init_voronoi_graph(vgraph_t *vgraph,
			       const nb_mesh_t *const mesh)
{
	vgraph->N = vcn_mesh_get_N_vtx(mesh);
	uint32_t N_adj = vcn_mesh_get_N_edg(mesh);
	
	uint32_t size1 = vgraph->N * sizeof(*(vgraph->N_adj));
	uint32_t size2 = vgraph->N * sizeof(*(vgraph->adj));
	uint32_t memsize = size1 + size2 + 
		2 * N_edg * sizeof(**(vgraph->adj));
	char *memblock = malloc(memsize);

	vgraph->N_adj = (void*) memblock;
	vgraph->adj = (void*) (memblock +  size1);

	count_vgraph_adj(vgraph, mesh);
	set_vgraph_adj_mem(vgraph, memblock + size1 + size2);
	set_vgraph_adj(vgraph, mesh);
}

static void count_vgraph_adj(vgraph_t *vgraph, const nb_mesh_t *const mesh)
{
	/* AQUI VOY */
	memset(vgraph->N_adj, 0,  vgraph->N * sizeof(*(vgraph->N_adj)));

	uint16_t iter_size = nb_iterator_get_memsize();
	nb_iterator_t *iter = alloca(iter_size);
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, mesh->ht_edge);
	while (nb_iterator_has_more(iter)) {
		msh_edge_t *edge = (msh_edge_t*) nb_iterator_get_next(iter);
		uint32_t id1 = *(uint32_t*)((void**)edge->v1->attr)[0];
		uint32_t id2 = *(uint32_t*)((void**)edge->v2->attr)[0];
		vgraph->N_adj[id1] += 1;
		vgraph->N_adj[id2] += 1;
	}
	nb_iterator_finish(iter);
}

static void set_vgraph_adj_mem(vgraph_t *vgraph, char *memblock)
{
	for (uint32_t i = 0; i < vgraph->N; i++) {
		vgraph->adj[i] = (void*) memblock;
		memblock += vgraph->N_adj[i] * sizeof(**(vgraph->adj));
	}
}

static void set_vgraph_adj(vgraph_t *vgraph, const nb_mesh_t *const mesh)
{
	memset(vgraph->N_adj, 0,  vgraph->N * sizeof(*(vgraph->N_adj)));

	uint16_t iter_size = nb_iterator_get_memsize();
	nb_iterator_t *iter = alloca(iter_size);
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, mesh->ht_edge);
	while (nb_iterator_has_more(iter)) {
		msh_edge_t *edge = (msh_edge_t*) nb_iterator_get_next(iter);
		bool without_cc = !adj_is_cocircular(edge);
		if (without_cc) {
			uint32_t id1 = *(uint32_t*)((void**)edge->v1->attr)[0];
			uint32_t id2 = *(uint32_t*)((void**)edge->v2->attr)[0];

			vgraph->adj[id1][vgraph->N_adj[id1]] = edge;
			vgraph->adj[id2][vgraph->N_adj[id2]] = edge;

			vgraph->N_adj[id1] += 1;
			vgraph->N_adj[id2] += 1;			
		}
	}
	nb_iterator_finish(iter);
}

static bool adj_is_cocircular(const msh_edge_t *const edg)
{
	bool out;
	if (medge_is_boundary(edg)) {
		out = edge_has_trg_with_cl_ccenter(edg);
	} else {
		msh_vtx_t *v4 = medge_get_opposite_vertex(edg->t2, edg);
		out = nb_utils2D_pnt_is_cocircular(edg->t1->v1->x,
						   edg->t1->v2->x,
						   edg->t1->v3->x,
						   v4->x);
	}
	return out;
}

static bool edge_has_trg_with_cl_ccenter(const msh_edge_t *const edg)
/* Has a colinear circumcenter of the adjacent triangles? */
{
	bool is_cl;
	if (NULL != edg->t1)
		is_cl = edge_is_encroached(edg, edg->t1);
	else
		is_cl = false;

	if (!is_cl) {
		if (NULL != edg->t2)
			is_cl = edge_is_encroached(edg, edg->t2);
	}
	return is_cl;
}

static bool edge_is_encroached(const msh_edge_t *edg, const msh_trg_t *trg)
{
	msh_vtx_t  *p;
	if (trg->s1 == edg)
		p = trg->v3;
	else if (trg->s2 == edg)
		p = trg->v1;
	else
		p = trg->v2;
	return vcn_utils2D_pnt_lies_in_diametral_circle(edg->v1->x,
							edg->v2->x,
							p->x);
}

static void finish_voronoi_graph(vgraph_t *vgraph)
{
	free(vgraph->N_adj);
}

static void set_voronoi(nb_mshpoly_t *poly,
			const vgraph_t *const vgraph)
{
	count_voronoi_data(poly, vgraph);

	poly->N_vtx = mesh->N_input_vtx;
	poly->N_sgm = mesh->N_input_sgm;

	set_arrays_memory(poly);

	set_nodes(poly, mesh, cc_map);
	set_edges(poly, mesh);

	uint32_t memsize = set_N_sides(poly, mesh);
	set_mem_of_adj_and_ngb(poly, memsize);
	set_adj_and_ngb(poly, mesh);

	set_elem_vtx(poly, mesh);

	memsize = set_N_nod_x_sgm(poly, mesh);
	set_mem_of_nod_x_sgm(poly, memsize);
	set_nod_x_sgm(poly, mesh);
}

static void count_voronoi_data(nb_mshpoly_t *poly,
			       const vgraph_t *const vgraph)
{
	uint32_t N_duplicate_edg = 0;
	uint32_t N_boundary_cells = 0;
	for (uint32_t i = 0; i < vgraph->N; i++) {
		N_duplicate_edg += vgraph->N_adj[i];
		if (cell_is_boundary(vgraph, i))
			N_boundary_cells += 1;
	}
	poly->N_nod = N_duplicate_edg / 2;
	poly->N_edg = N_duplicate_edg / 2 + 2 * N_boundary_cells;
	poly->N_elems = vgraph->N;
}

static bool cell_is_boundary(const vgraph_t *const vgraph,
			     uint32_t i)
{
/* NOT ALWAYS WORK (Colineal boundary ccenters */
	bool out = false;
	for (uint16_t j = 0; j < vgraph->N_adj[i]; j++) {
		if (medge_is_boundary(vgraph->adj[i][j])) {
			out = true;
			break;
		}
	}
	return out;
}

static void set_nodes(nb_mshpoly_t *poly,
		      const nb_mesh_t *const mesh,
		      const uint32_t *cc_map);
static void set_edges(nb_mshpoly_t *poly,
		      const nb_mesh_t *const mesh);
static uint32_t set_N_sides(nb_mshpoly_t *poly,
			   const nb_mesh_t *const mesh);
static void set_adj_and_ngb(nb_mshpoly_t *poly,
			    const nb_mesh_t *const mesh);
static void set_elem_vtx(nb_mshpoly_t *poly,
			 const nb_mesh_t *const mesh);
static uint32_t set_N_nod_x_sgm(nb_mshpoly_t *poly,
				const nb_mesh_t *const mesh);
static void set_nod_x_sgm(nb_mshpoly_t *poly,
			  const nb_mesh_t *const mesh);

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
