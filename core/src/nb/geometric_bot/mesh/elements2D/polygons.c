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
	int8_t *type; /* 0: Interior, 1: Subsegment */
	uint16_t *N_adj;
	msh_edge_t ***adj;
} vgraph_t;

typedef struct {
	uint32_t N_trg_in;  /* # interior trg */
	uint32_t N_vtx_in;  /* # interior vtx */
	uint32_t N_vtx_out; /* # vtx on input sgm  */
	uint32_t N_edg_in;  /* # interior edges */
	uint32_t N_edg_out; /* # edg on input sgm */
	uint32_t N_cc_in;   /* # Interior edg joining
			     *   cocircular interior trg
			     */
	uint32_t *vtx_map;
	uint32_t *trg_map;
} vinfo_t;

static void set_arrays_memory(nb_mshpoly_t *poly);
static void copy_nodes(nb_mshpoly_t* poly,
		       const nb_mshpoly_t *const src_poly);
static void copy_edges(nb_mshpoly_t* poly,
		       const nb_mshpoly_t *const src_poly);
static void copy_centroids(nb_mshpoly_t* poly,
			   const nb_mshpoly_t *const src_poly);
static void copy_N_adj(nb_mshpoly_t* poly,
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

static void init_voronoi_info(vinfo_t *vinfo,
			      const nb_mesh_t *const mesh);
static void init_trg_cc_map(uint32_t *trg_cc_map, uint32_t Nt);
static void init_voronoi_graph(vgraph_t *vgraph, vinfo_t *vinfo,
			       uint32_t *trg_cc_map,
			       const nb_mesh_t *const mesh);
static void create_mapping(vinfo_t *vinfo,
			   const vgraph_t *const vgraph,
			   const nb_mesh_t *const mesh,
			   const uint32_t *trg_cc_map);
static void count_vgraph_adj(vgraph_t *vgraph, vinfo_t *vinfo,
			     const nb_mesh_t *const mesh);
static void set_vgraph_adj_mem(vgraph_t *vgraph, char *memblock);
static void set_vgraph_adj(vgraph_t *vgraph, vinfo_t *vinfo,
			   uint32_t *trg_cc_map,
			   const nb_mesh_t *const mesh);
static bool adj_is_cocircular(const msh_edge_t *const edg);
static void update_cc_map(const msh_edge_t *edge, bool is_cc,
			  uint32_t *trg_cc_map, uint32_t N_trg);
static void insert_edg_as_adj(vgraph_t *vgraph, const msh_edge_t *edge,
			      bool is_cc);
static void insert_adj_sorted_by_angle(vgraph_t *vgraph, uint16_t igraph,
				       msh_edge_t *edge);
static void counting_edg_in_vinfo(vinfo_t *vinfo, const vgraph_t *vgraph,
				  const msh_edge_t *edge, bool is_cc);
static void finish_voronoi_info(vinfo_t *vinfo);
static void finish_voronoi_graph(vgraph_t *vgraph);
static void set_voronoi(nb_mshpoly_t *poly,
			const vgraph_t *const vgraph,
			const vinfo_t *const vinfo,
			const nb_mesh_t *const mesh);
static void set_quantities(nb_mshpoly_t *poly,
			   const vinfo_t *const vinfo,
			   const nb_mesh_t *const mesh);
static void set_nodes_and_centroids(nb_mshpoly_t *poly,
				    const vgraph_t *const vgraph,
				    const vinfo_t *const vinfo,
				    const nb_mesh_t *const mesh);
static void set_edges(nb_mshpoly_t *poly,
		      const vgraph_t *const vgraph,
		      const vinfo_t *const vinfo,
		      const nb_mesh_t *const mesh);

static void process_interior_edge(nb_mshpoly_t *poly,
				  const vgraph_t *const vgraph,
				  const vinfo_t *const vinfo,
				  const msh_edge_t *const edg,
				  uint32_t iedge);
static void set_N_adj(nb_mshpoly_t *poly,
		      const vgraph_t *const vgraph,
		      const vinfo_t *const vinfo);

static void set_adj_and_ngb(nb_mshpoly_t *poly,
			    const vgraph_t *const vgraph,
			    const vinfo_t *const vinfo);
static uint16_t add_adj_and_ngb(nb_mshpoly_t *poly,
				const vgraph_t *const vgraph,
				const vinfo_t *const vinfo,
				uint32_t i, uint16_t j, uint16_t id_adj);
static msh_vtx_t *get_partner(const vgraph_t *const vgraph,
			      const vinfo_t *const vinfo,
			      uint32_t i, uint16_t j);
static msh_trg_t *get_prev_trg(const vgraph_t *const vgraph,
			       const vinfo_t *const vinfo,
			       uint32_t i, uint16_t j);

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
		copy_centroids(poly, src_poly);
		copy_N_adj(poly, src_poly);

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
	uint32_t cen_size = poly->N_elems * sizeof(*(poly->cen));
	uint32_t N_adj_size = poly->N_elems * sizeof(*(poly->N_adj));
	uint32_t adj_size = poly->N_elems * sizeof(*(poly->adj));
	uint32_t ngb_size = poly->N_elems * sizeof(*(poly->ngb));
	uint32_t elem_vtx_size = poly->N_vtx * sizeof(*(poly->elem_vtx));
	uint32_t N_nod_x_sgm_size = poly->N_sgm * sizeof(*(poly->N_nod_x_sgm));
	uint32_t nod_x_sgm_size = poly->N_sgm * sizeof(*(poly->nod_x_sgm));

	uint32_t size = nod_size + edg_size + cen_size + N_adj_size +
		adj_size + ngb_size + elem_vtx_size +
		N_nod_x_sgm_size + nod_x_sgm_size;

	char *memblock = malloc(size);

	poly->nod = (void*) memblock;
	poly->edg = (void*) ((char*)(poly->nod) + nod_size);
	poly->cen = (void*) ((char*)(poly->edg) + edg_size);
	poly->N_adj = (void*) ((char*)(poly->cen) + cen_size);
	poly->adj = (void*) ((char*)(poly->N_adj) + N_adj_size);
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

static void copy_centroids(nb_mshpoly_t* poly,
			   const nb_mshpoly_t *const src_poly)
{
	memcpy(poly->cen, src_poly->cen,
	       2 * poly->N_elems * sizeof(*(poly->cen)));
}

static void copy_N_adj(nb_mshpoly_t* poly, const nb_mshpoly_t *const src_poly)
{
	memcpy(poly->N_adj, src_poly->N_adj,
	       poly->N_elems * sizeof(*(poly->N_adj)));
}

static uint32_t get_size_of_adj_and_ngb(const nb_mshpoly_t *const poly)
{
	uint32_t size = 0;
	for (uint32_t i = 0; i < poly->N_sgm; i++)
		size += poly->N_adj[i] * sizeof(**(poly->adj));
	return 2 * size;
}

static void set_mem_of_adj_and_ngb(nb_mshpoly_t *poly, uint32_t memsize)
{
	char *memblock1 = malloc(memsize);
	char *memblock2 = memblock1 + memsize / 2;
	for (uint32_t i = 0; i < poly->N_sgm; i++) {
		poly->adj[i] = (void*) memblock1;
		poly->ngb[i] = (void*) memblock2;
		memblock1 += poly->N_adj[i] * sizeof(**(poly->adj));
		memblock2 += poly->N_adj[i] * sizeof(**(poly->ngb));
	}
}

static void copy_adj_and_ngb(nb_mshpoly_t* poly,
			     const nb_mshpoly_t *const src_poly)
{
	for (int i = 0; i < poly->N_sgm; i++) {
		memcpy(&(poly->adj[i]), &(src_poly->adj[i]),
		       poly->N_adj[i] * sizeof(**(poly->adj)));
		memcpy(&(poly->ngb[i]), &(src_poly->ngb[i]),
		       poly->N_adj[i] * sizeof(**(poly->ngb)));
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

void nb_mshpoly_load_from_mesh(nb_mshpoly_t *mshpoly, nb_mesh_t *mesh)
{
	if (vcn_mesh_get_N_trg(mesh) > 0) {
		/* POR IMPLEMENTAR */mesh_split_trg_formed_by_input_sgm(mesh);

		mesh_alloc_vtx_ids((vcn_mesh_t*)mesh);
		mesh_alloc_trg_ids((vcn_mesh_t*)mesh);
	
		vinfo_t vinfo;
		init_voronoi_info(&vinfo, mesh);

		uint32_t *trg_cc_map = malloc(vcn_mesh_get_N_trg(mesh) * sizeof(uint32_t));
		init_trg_cc_map(trg_cc_map);

		vgraph_t vgraph;
		init_voronoi_graph(&vgraph, &vinfo, trg_cc_map, mesh);

		create_mapping(&vinfo, &vgraph, mesh, trg_cc_map);
		free(trg_cc_map);

		set_voronoi(mshpoly, &vgraph, &vinfo, mesh);

		finish_voronoi_info(&vinfo);
		finish_voronoi_graph(&vgraph);

		mesh_free_vtx_ids((vcn_mesh_t*)mesh);
		mesh_free_trg_ids((vcn_mesh_t*)mesh);
	}
}

static void init_voronoi_info(vinfo_t *vinfo,
			      const nb_mesh_t *const mesh)
{
	memset(vinfo, 0,  sizeof(*vinfo));
	uint32_t Nv = vcn_mesh_get_N_vtx();
	uint32_t Nt = vcn_mesh_get_N_trg();
	uint32_t size1 = Nv * sizeof(uint32_t);
	uint32_t size2 = Nt * sizeof(uint32_t);
	char *memblock = malloc(2 * size1 + size2);
	vinfo->vtx_map = (void*) memblock;
	vinfo->trg_map = (void*) (memblock + size1);

	for (uint32_t i = 0; i < Nv; i++)
		vinfo->vtx_map[i] = Nv;

	for (uint32_t i = 0; i < Nt; i++)
		vinfo->trg_map[i] = Nt;
}

static void init_trg_cc_map(uint32_t *trg_cc_map, uint32_t Nt)
{
	for (uint32_t i = 0; i < Nt; i++)
		trg_cc_map[i] = i;
}

static void init_voronoi_graph(vgraph_t *vgraph, vinfo_t *vinfo,
			       uint32_t *trg_cc_map,
			       const nb_mesh_t *const mesh)
{
	vgraph->N = vcn_mesh_get_N_vtx(mesh);
	
	uint32_t size1 = vgraph->N * sizeof(*(vgraph->type));
	uint32_t size2 = vgraph->N * sizeof(*(vgraph->N_adj));
	uint32_t size3 = vgraph->N * sizeof(*(vgraph->adj));
	uint32_t memsize = size1 + size2 + size3 +
		2 * N_edg * sizeof(**(vgraph->adj));
	char *memblock = malloc(memsize);

	vgraph->type = (void*) memblock;
	vgraph->N_adj = (void*) (memblock +  size1);
	vgraph->adj = (void*) (memblock +  size1 + size2);

	count_vgraph_adj(vgraph, vinfo, mesh);
	set_vgraph_adj_mem(vgraph, memblock + size1 + size2 + size3);
	set_vgraph_adj(vgraph, vinfo, trg_cc_map, mesh);
}

static void create_mapping(vinfo_t *vinfo,
			   const vgraph_t *const vgraph,
			   const nb_mesh_t *const mesh,
			   const uint32_t *trg_cc_map)
{
	uint32_t ielem = 0;
	uint32_t inode = 0;

	vcn_bins2D_iter_t* biter = vcn_bins2D_iter_create();
	vcn_bins2D_iter_set_bins(biter, mesh->ug_vtx);
	while (vcn_bins2D_iter_has_more(biter)) {
		const msh_vtx_t* vtx = vcn_bins2D_iter_get_next(biter);
		uint32_t id = *(uint32_t*)((void**)vtx->attr)[0];
		bool vtx_is_interior = 0 == vgraph->type[id];
		if (vtx_is_interior) {
			vinfo->vtx_map[id] = ielem;
			ielem += 1;
		} else {
			vinfo->vtx_map[id] = inode;
			inode += 1;
		}
	}
	vcn_bins2D_iter_destroy(iter);

	uint16_t iter_size = nb_iterator_get_memsize();
	nb_iterator_t *iter = alloca(iter_size);
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, mesh->ht_trg);
	while (nb_iterator_has_more(iter)) {
		const msh_trg_t *trg = nb_iterator_get_next(iter);
		uint32_t id = *(uint32_t*)((void**)trg->attr)[0];
		uint32_t v1 = *(uint32_t*)((void**)trg->v1->attr)[0];
		uint32_t v2 = *(uint32_t*)((void**)trg->v2->attr)[0];
		uint32_t v3 = *(uint32_t*)((void**)trg->v3->attr)[0];
		bool trg_is_interior = (0 == vgraph->type[v1]) &&
			(0 == vgraph->type[v2]) && (0 == vgraph->type[v3]);
		if (trg_is_interior) {
			if (id != trg_cc_map[id]) {
				cc_id = trg_cc_map[id];
				vinfo->trg_map[id] = vinfo->trg_map[cc_id];
			} else {
				vinfo->trg_map[id] = inode;
				inode += 1;
			}
		}
	}
	nb_iterator_finish(iter);
	
}

static void count_vgraph_adj(vgraph_t *vgraph, vinfo_t *vinfo,
			     const nb_mesh_t *const mesh)
{
	memset(vgraph->type, 0,  vgraph->N * sizeof(*(vgraph->type)));
	memset(vgraph->N_adj, 0,  vgraph->N * sizeof(*(vgraph->N_adj)));

	uint16_t iter_size = nb_iterator_get_memsize();
	nb_iterator_t *iter = alloca(iter_size);
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, mesh->ht_edge);
	while (nb_iterator_has_more(iter)) {
		const msh_edge_t *edge = nb_iterator_get_next(iter);
		uint32_t id1 = *(uint32_t*)((void**)edge->v1->attr)[0];
		uint32_t id2 = *(uint32_t*)((void**)edge->v2->attr)[0];

		vgraph->N_adj[id1] += 1;
		vgraph->N_adj[id2] += 1;

		if (medge_is_subsgm(edge)) {
			if (0 == vgraph->type[id1]) {
				vgraph->type[id1] = 1;
				vinfo->N_vtx_out += 1;
			}
			if (0 == vgraph->type[id2]) {
				vgraph->type[id2] = 1;
				vinfo->N_vtx_out += 1;
			}
		}
	}
	nb_iterator_finish(iter);

	vinfo->N_vtx_in = vgraph->N - vinfo->N_vtx_out;
}

static void set_vgraph_adj_mem(vgraph_t *vgraph, char *memblock)
{
	for (uint32_t i = 0; i < vgraph->N; i++) {
		vgraph->adj[i] = (void*) memblock;
		memblock += vgraph->N_adj[i] * sizeof(**(vgraph->adj));
	}
}

static void set_vgraph_adj(vgraph_t *vgraph, vinfo_t *vinfo,
			   uint32_t *trg_cc_map,
			   const nb_mesh_t *const mesh)
{
	memset(vgraph->N_adj, 0,  vgraph->N * sizeof(*(vgraph->N_adj)));

	uint16_t iter_size = nb_iterator_get_memsize();
	nb_iterator_t *iter = alloca(iter_size);
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, mesh->ht_edge);
	while (nb_iterator_has_more(iter)) {
		const msh_edge_t *edge = nb_iterator_get_next(iter);

		bool cc = adj_is_cocircular(edge);
		update_cc_map(edge, cc, trg_cc_map,
			      vcn_mesh_get_N_trg(mesh));
		insert_edg_as_adj(vgraph, edge, cc);
		counting_edg_in_vinfo(vinfo, vgraph, edge, cc);
	}
	nb_iterator_finish(iter);
	
	vinfo->N_trg_in /= 3;
}

static bool adj_is_cocircular(const msh_edge_t *const edg)
{
	bool out;
	if (medge_is_boundary(edg)) {
		out = false;
	} else {
		msh_vtx_t *v4 = medge_get_opposite_vertex(edg->t2, edg);
		out = nb_utils2D_pnt_is_cocircular(edg->t1->v1->x,
						   edg->t1->v2->x,
						   edg->t1->v3->x,
						   v4->x);
	}
	return out;
}

static void update_cc_map(const msh_edge_t *edge, bool is_cc,
			  uint32_t *trg_cc_map, uint32_t N_trg)
{
	if (is_cc) {
		uint32_t id1 = *(uint32_t*)((void**)edge->t1->attr)[0];
		uint32_t id2 = *(uint32_t*)((void**)edge->t2->attr)[0];
		if (id1 == trg_cc_map[id1] && id2 == trg_cc_map[id2]) {
			trg_cc_map[id1] = id2;
		} else if (id1 == trg_cc_map[id1] && id2 != trg_cc_map[id2]) {
			trg_cc_map[id1] = trg_cc_map[id2];
		} else if (id1 != trg_cc_map[id1] && id2 == trg_cc_map[id2]) {
			trg_cc_map[id2] = trg_cc_map[id1];
		} else {
			uint32_t new_cc = trg_cc_map[id1];
			uint32_t old_cc = trg_cc_map[id2];
			for (uint32_t i = 0; i < N_trg; i++) {
				if (old_cc == trg_cc_map[i])
					trg_cc_map[i] = new_cc;
			}				
		}
	}
}

static void insert_edg_as_adj(vgraph_t *vgraph, const msh_edge_t *edge,
			      bool is_cc)
{
	uint32_t id1 = *(uint32_t*)((void**)edge->v1->attr)[0];
	uint32_t id2 = *(uint32_t*)((void**)edge->v2->attr)[0];
	bool not_interior_cocircular_edg = !is_cc ||
		1 == vgraph->type[id1] || 1 == vgraph->type[id2];
	bool semi_interior =
		0 == vgraph->type[id1] || 0 == vgraph->type[id2];
	if (not_interior_cocircular_edg && semi_interior) {
		insert_adj_sorted_by_angle(vgraph, id1, edge);
		insert_adj_sorted_by_angle(vgraph, id2, edge);
	}
}

static void insert_adj_sorted_by_angle(vgraph_t *vgraph, uint16_t igraph,
				       msh_edge_t *edge)
{
	msh_vtx_t *v1, v2;
	if (igraph == *(uint32_t*)((void**)edge->v1->attr)[0]) {
		v1 = edge->v1;
		v2 = edge->v2;
	} else {
		v1 = edge->v2;
		v2 = edge->v1;
	}

	uint32_t id = vgraph->N_adj[igraph];
	double angle_id;
	if (0 < id) {
		double x = v2->x[0] - v1->x[0];
		double y = v2->x[1] - v1->x[1];
		angle_id = atan2(y, x);
	}

	uint16_t j = 0;
	while (j < id) {
		msh_vtx_t *v2 = medge_get_partner_vtx(vgraph->adj[igraph][j], v1);
		double x = v2->x[0] - v1->x[0];
		double y = v2->x[1] - v1->x[1];
		double angle_j = atan2(y, x);
		if (angle_j > angle_id) {
			msh_edge_t *aux = vgraph->adj[igraph][j];
			vgraph->adj[igraph][j] = edge;
			edge = aux;
			angle_id = angle_j;
		}
	}
	vgraph->adj[igraph][id] = edge;
	vgraph->N_adj[igraph] += 1;
}

static void counting_edg_in_vinfo(vinfo_t *vinfo, const vgraph_t *vgraph,
				  const msh_edge_t *edge, bool is_cc)
{
	uint32_t id1 = *(uint32_t*)((void**)edge->v1->attr)[0];
	uint32_t id2 = *(uint32_t*)((void**)edge->v2->attr)[0];
	if (1 == vgraph->type[id1] && 1 == vgraph->type[id2]) {
		vinfo->N_edg_out += 1;
	} else if (0 == vgraph->type[id1] && 0 == vgraph->type[id2]) {
		vinfo->N_edg_in += 1;

		msh_vtx_t *opp_t1 = mtrg_get_opposite_vertex(edge->t1, edge);
		msh_vtx_t *opp_t2 = mtrg_get_opposite_vertex(edge->t2, edge);
		uint32_t opp1 = *(uint32_t*)((void**)opp_t1->attr)[0];
		uint32_t opp2 = *(uint32_t*)((void**)opp_t2->attr)[0];
		bool trg_are_interior =
			0 == vgraph->type[opp1] && 0 == vgraph->type[opp2];
		if (cc && trg_are_interior)
			vinfo->N_cc_in += 1;

		if (0 == vgraph->type[opp1])
			vinfo->N_trg_in += 1;

		if (0 == vgraph->type[opp2])
			vinfo->N_trg_in += 1;
	}
}

static void finish_voronoi_info(vinfo_t *vinfo)
{
	free(vinfo->vtx_map);
}

static void finish_voronoi_graph(vgraph_t *vgraph)
{
	free(vgraph->type);
}

static void set_voronoi(nb_mshpoly_t *poly,
			const vgraph_t *const vgraph,
			const vinfo_t *const vinfo,
			const nb_mesh_t *const mesh)
{
	set_quantities(poly, &vinfo, mesh);
	
	set_arrays_memory(poly);

	set_nodes_and_centroids(poly, vgraph, vinfo, mesh);
	set_edges(poly, vgraph, vinfo, mesh);
	set_N_adj(poly, vgraph, vinfo);
	
	uint32_t memsize = get_size_of_adj_and_ngb(poly);
	set_mem_of_adj_and_ngb(poly, memsize);
	set_adj_and_ngb(poly, vgraph, vinfo);

	/* POR IMPLEMENTAR */set_elem_vtx(poly);
	/* POR IMPLEMENTAR */set_N_nod_x_sgm(poly);

	memsize = get_size_of_nod_x_sgm(poly);
	set_mem_of_nod_x_sgm(poly, memsize);
	/* POR IMPLEMENTAR */set_nod_x_sgm(poly);
}


static void set_quantities(nb_mshpoly_t *poly,
			   const vinfo_t *const vinfo,
			   const nb_mesh_t *const mesh)
{
	poly->N_nod = vinfo->N_trg_in + vinfo->N_vtx_out - vinfo->N_cc_in;
	poly->N_edg = vinfo->N_edg_int + vinfo->N_edg_out - vinfo->N_cc_in;
	poly->N_elems = vinfo->N_vtx_in;

	poly->N_vtx = mesh->N_input_vtx;
	poly->N_sgm = mesh->N_input_sgm;
}

static void set_nodes_and_centroids(nb_mshpoly_t *poly,
				    const vgraph_t *const vgraph,
				    const vinfo_t *const vinfo,
				    const nb_mesh_t *const mesh)
{
	vcn_bins2D_iter_t* biter = vcn_bins2D_iter_create();
	vcn_bins2D_iter_set_bins(biter, mesh->ug_vtx);
	while (vcn_bins2D_iter_has_more(biter)) {
		const msh_vtx_t* vtx = vcn_bins2D_iter_get_next(biter);
		uint32_t id = *(uint32_t*)((void**)vtx->attr)[0];
		bool vtx_is_interior = 0 == vgraph->type[id];
		if (vtx_is_interior) {
			uint32_t ielem = vinfo->vtx_map[id];
			memcpy(&(poly->cen[ielem * 2]), vtx->x,
			       2 * sizeof(*(vtx->x)));
		} else {
			uint32_t inode = vinfo->vtx_map[id];
			memcpy(&(poly->nod[inode * 2]), vtx->x,
			       2 * sizeof(*(vtx->x)));
		}
	}
	vcn_bins2D_iter_destroy(iter);

	uint16_t iter_size = nb_iterator_get_memsize();
	nb_iterator_t *iter = alloca(iter_size);
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, mesh->ht_trg);
	while (nb_iterator_has_more(iter)) {
		const msh_trg_t *trg = nb_iterator_get_next(iter);
		uint32_t id = *(uint32_t*)((void**)trg->attr)[0];
		uint32_t v1 = *(uint32_t*)((void**)trg->v1->attr)[0];
		uint32_t v2 = *(uint32_t*)((void**)trg->v2->attr)[0];
		uint32_t v3 = *(uint32_t*)((void**)trg->v3->attr)[0];
		bool trg_is_interior = (0 == vgraph->type[v1]) &&
			(0 == vgraph->type[v2]) && (0 == vgraph->type[v3]);
		if (trg_is_interior) {
			uint32_t inode = vinfo->trg_map[id];
			double circumcenter[2];
			vcn_utils2D_get_circumcenter(trg->v1->x,
						     trg->v2->x,
						     trg->v3->x,
						     circumcenter);
			memcpy(&(poly->nod[inode * 2]), circumcenter,
			       2 * sizeof(*(trg->v1->x)));
		}
	}
	nb_iterator_finish(iter);
}

static void set_edges(nb_mshpoly_t *poly,
		      const vgraph_t *const vgraph,
		      const vinfo_t *const vinfo,
		      const nb_mesh_t *const mesh)
{
	uint32_t iedge = 0;

	uint16_t iter_size = nb_iterator_get_memsize();
	nb_iterator_t *iter = alloca(iter_size);
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, mesh->ht_edge);
	while (nb_iterator_has_more(iter)) {
		const msh_edge_t *edg = nb_iterator_get_next(iter);
		uint32_t v1 = *(uint32_t*)((void**)edg->v1->attr)[0];
		uint32_t v2 = *(uint32_t*)((void**)edg->v2->attr)[0];
		if (0 == vgraph->type[v1] && 0 == vgraph->type[v2]) {
			iedge = process_interior_edge(poly, vgraph, vinfo, edg);
			iedge += 1;
		} else if (1 == vgraph->type[v1] && 1 == vgraph->type[v2]) {
			/* Process not interior edges */
			poly->edg[iedge * 2] = vinfo->vtx_map[v1];
			poly->edg[iedge*2+1] = vinfo->vtx_map[v2];
			iedge += 1;
		}
	}
	nb_iterator_finish(iter);
}

static void process_interior_edge(nb_mshpoly_t *poly,
				  const vgraph_t *const vgraph,
				  const vinfo_t *const vinfo,
				  const msh_edge_t *const edg,
				  uint32_t iedge)
{
	msh_vtx_t *opp_t1 = mtrg_get_opposite_vertex(edg->t1, edg);
	msh_vtx_t *opp_t2 = mtrg_get_opposite_vertex(edg->t2, edg);
	uint32_t opp1 = *(uint32_t*)((void**)opp_t1->attr)[0];
	uint32_t opp2 = *(uint32_t*)((void**)opp_t2->attr)[0];
	bool trg_are_interior =
		0 == vgraph->type[opp1] && 0 == vgraph->type[opp2];

	uint32_t t1 = *(uint32_t*)((void**)edg->t1->attr)[0];
	uint32_t t2 = *(uint32_t*)((void**)edg->t2->attr)[0];
	bool is_not_cc = vinfo->trg_map[t1] != vinfo->trg_map[t2];
	if (is_not_cc && trg_are_interior(edg, vgraph)) {
		/* Interior triangles (No cocircular) */
		poly->edg[iedge * 2] = vinfo->trg_map[t1];
		poly->edg[iedge*2+1] = vinfo->trg_map[t2];		
	} else {
		/* At least a triangle is not interior */
		uint32_t node1;
		if (0 == vgraph->type[opp1])
			node1 = vinfo->trg_map[t1];
		else
			node1 = vinfo->vtx_map[opp1];
		uint32_t node2;
		if (0 == vgraph->type[opp2])
			node2 = vinfo->trg_map[t2];
		else
			node2 = vinfo->vtx_map[opp2];

		poly->edg[iedge * 2] = node1;
		poly->edg[iedge*2+1] = node2;			
	}
}

static void set_N_adj(nb_mshpoly_t *poly,
		      const vgraph_t *const vgraph,
		      const vinfo_t *const vinfo)
{
	for (uint32_t i = 0; i < vgraph->N; i++) {
		if (0 == vgraph->type[i]) {
			uint32_t elem_id = vinfo->vtx_map[i];
			poly->N_adj[elem_id] = vgraph->N_adj[i];
		}
	}
}

static void set_adj_and_ngb(nb_mshpoly_t *poly,
			    const vgraph_t *const vgraph,
			    const vinfo_t *const vinfo)
{
	for (uint32_t i = 0; i < vgraph->N; i++) {
		if (0 == vgraph->type[i]) {
			uint32_t elem_id = vinfo->vtx_map[i];
			uint16_t id_adj = 0;
			for (uint16_t j = 0; j < vgraph->N_adj[i]; j++) {
				id_adj = add_adj_and_ngb(poly, vgraph,
							 vinfo, i, j,
							 id_adj);
			}
		}
	}
}

static uint16_t add_adj_and_ngb(nb_mshpoly_t *poly,
				const vgraph_t *const vgraph,
				const vinfo_t *const vinfo,
				uint32_t i, uint16_t j, uint16_t id_adj)
{
	uint16_t id_prev = (j + vgraph->N_adj[i] - 1) % vgraph->N_adj[i];
	msh_vtx_t *v1 = get_partner(vgraph, i, id_prev);
	msh_vtx_t *v2 = get_partner(vgraph, i, j);

	uint32_t id1 = *(uint32_t*)((void**)v1->attr)[0];
	uint32_t id2 = *(uint32_t*)((void**)v2->attr)[0];
	if (0 == vgraph->type[id2]) {
		if (0 == vgraph->type[id1]) {
			/* Interior trg case */
			uint32_t trg_id = get_prev_trg(vgraph, i, j);
			poly->adj[i][id_adj] = vinfo->trg_map[trg_id];
		} else {
			/* A node in the boundary */
			poly->adj[i][id_adj] = vinfo->vtx_map[id1];
		}
		poly->ngb[i][id_adj] = v2;
	} else {
		poly->adj[i][id_adj] = vinfo->vtx_map[id_ngb];
		poly->ngb[i][id_adj] = poly->N_elems;
		id_adj += 1;
	}
}

static msh_vtx_t *get_partner(const vgraph_t *const vgraph,
			      const vinfo_t *const vinfo,
			      uint32_t i, uint16_t j)
{
	msh_vtx_t *out = NULL;
	if (j < vgraph->N_adj[i]) {
		msh_edge_t *edge = vgraph->adj[i][j];
		if (NULL != edge) {
			if (i == *(uint32_t*)((void**)edge->v1->attr)[0])
				out = edge->v2;
			else
				out = edge->v1;
		}
	}
	return out;
}

static msh_trg_t *get_prev_trg(const vgraph_t *const vgraph,
			       const vinfo_t *const vinfo,
			       uint32_t i, uint16_t j)
{
/* AQUI VOY */
}

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
