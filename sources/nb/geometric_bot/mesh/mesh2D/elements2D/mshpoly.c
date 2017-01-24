#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "nb/memory_bot.h"
#include "nb/geometric_bot.h"

#include "../../tessellator2D_structs.h"
#include "../../ruppert/ruppert.h"
#include "mshpoly_struct.h"


#define POW2(a) ((a)*(a))

typedef struct {
	uint32_t N;
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

static void copy_nodes(nb_mshpoly_t* poly,
		       const nb_mshpoly_t *const src_poly);
static void copy_edges(nb_mshpoly_t* poly,
		       const nb_mshpoly_t *const src_poly);
static void copy_centroids(nb_mshpoly_t* poly,
			   const nb_mshpoly_t *const src_poly);
static void copy_N_adj(nb_mshpoly_t* poly,
			const nb_mshpoly_t *const src_poly);
static uint32_t get_size_of_adj_and_ngb(const nb_mshpoly_t *const poly);
static void copy_adj_and_ngb(nb_mshpoly_t* poly,
			     const nb_mshpoly_t *const src_poly);

static void copy_elem_vtx(nb_mshpoly_t* poly,
			  const nb_mshpoly_t *const src_poly);
static void copy_N_nod_x_sgm(nb_mshpoly_t* poly,
			     const nb_mshpoly_t *const src_poly);
static uint32_t get_size_of_nod_x_sgm(const nb_mshpoly_t *const poly);
static void copy_nod_x_sgm(nb_mshpoly_t* poly,
			   const nb_mshpoly_t *const src_poly);
static void* nb_allocate_mem_poly(void);

static void init_voronoi_info(vinfo_t *vinfo,
			      const nb_tessellator2D_t *const mesh);
static void set_nodal_perm_to_nodes(nb_mshpoly_t *poly, const uint32_t *perm);
static void set_nodal_perm_to_edges(nb_mshpoly_t *poly, const uint32_t *perm);
static void set_nodal_perm_to_elems(nb_mshpoly_t *poly, const uint32_t *perm);
static void set_nodal_perm_to_invtx(nb_mshpoly_t *poly, const uint32_t *perm);
static void set_nodal_perm_to_insgm(nb_mshpoly_t *poly, const uint32_t *perm);
static void init_trg_cc_map(uint32_t *trg_cc_map, uint32_t Nt);
static void init_voronoi_graph(vgraph_t *vgraph, vinfo_t *vinfo,
			       uint32_t *trg_cc_map,
			       const nb_tessellator2D_t *const mesh);
static void create_mapping(vinfo_t *vinfo,
			   const vgraph_t *const vgraph,
			   const nb_tessellator2D_t *const mesh,
			   const uint32_t *trg_cc_map);
static bool trg_is_interior(const msh_trg_t *const trg,
			    const vgraph_t *const vgraph);
static void count_vgraph_adj(vgraph_t *vgraph,
			     const nb_tessellator2D_t *const mesh);
static void set_vgraph_adj_mem(vgraph_t *vgraph, char *memblock);
static void set_vgraph_adj(vgraph_t *vgraph, vinfo_t *vinfo,
			   uint32_t *trg_cc_map,
			   const nb_tessellator2D_t *const mesh);
static bool adj_is_cocircular(const msh_edge_t *const edg);
static void update_cc_map(const msh_edge_t *edge, bool is_cc,
			  uint32_t *trg_cc_map, uint32_t N_trg);
static void insert_edg_as_adj(vgraph_t *vgraph, const msh_edge_t *edge,
			      bool is_cc);
static void insert_adj_sorted_by_angle(vgraph_t *vgraph,
				       const msh_edge_t *edge,
				       const msh_vtx_t *v1,
				       const msh_vtx_t *v2);
static void counting_edg_in_vinfo(vinfo_t *vinfo, const vgraph_t *vgraph,
				  const msh_edge_t *edge, bool is_cc);
static void finish_voronoi_info(vinfo_t *vinfo);
static void finish_voronoi_graph(vgraph_t *vgraph);
static void set_voronoi(nb_mshpoly_t *poly,
			const vgraph_t *const vgraph,
			const vinfo_t *const vinfo,
			const nb_tessellator2D_t *const mesh);
static void set_quantities(nb_mshpoly_t *poly,
			   const vinfo_t *const vinfo,
			   const nb_tessellator2D_t *const mesh);
static void set_nodes_and_centroids(nb_mshpoly_t *poly,
				    const vgraph_t *const vgraph,
				    const vinfo_t *const vinfo,
				    const nb_tessellator2D_t *const mesh);
static void scale_vtx(double out[2], double in[2],
		      const nb_tessellator2D_t *mesh);
static void set_edges(nb_mshpoly_t *poly,
		      const vgraph_t *const vgraph,
		      const vinfo_t *const vinfo,
		      const nb_tessellator2D_t *const mesh);

static void process_interior_edge(nb_mshpoly_t *poly,
				  const vgraph_t *const vgraph,
				  const vinfo_t *const vinfo,
				  const msh_edge_t *const edg,
				  uint32_t iedge);
static void set_N_adj(nb_mshpoly_t *poly,
		      const vgraph_t *const vgraph,
		      const vinfo_t *const vinfo,
		      const nb_tessellator2D_t *const mesh);

static void set_adj_and_ngb(nb_mshpoly_t *poly,
			    const vgraph_t *const vgraph,
			    const vinfo_t *const vinfo,
			    const nb_tessellator2D_t *const mesh);
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
static void set_elem_vtx(nb_mshpoly_t *poly,
			 const vgraph_t *const vgraph,
			 const vinfo_t *const vinfo,
			 const nb_tessellator2D_t *const mesh);
static void set_N_nod_x_sgm(nb_mshpoly_t *poly,
			    const nb_tessellator2D_t *const mesh);
static void set_nod_x_sgm(nb_mshpoly_t *poly,
			  const vinfo_t *const vinfo,
			  const nb_tessellator2D_t *const mesh);
static void set_sgm_nodes(nb_mshpoly_t *poly,
			  const vinfo_t *const vinfo,
			  const nb_tessellator2D_t *const mesh,
			  uint32_t sgm_id);
static void assemble_sgm_wire(nb_mshpoly_t *poly,
			      const vinfo_t *const vinfo,
			      uint32_t sgm_id,
			      msh_edge_t *sgm_prev, msh_edge_t *sgm);
static void split_exterior_trg(nb_tessellator2D_t *mesh);
static void initialize_exterior_trg(const nb_tessellator2D_t *mesh,
				    nb_container_t *exterior_trg);
static bool is_exterior(const msh_trg_t *trg);
static bool have_all_nodes_in_sgm(const msh_trg_t *trg);
static void delete_exterior_trg(nb_tessellator2D_t *mesh,
				nb_container_t *exterior_trg);
static uint32_t get_nodes_N_adj(nb_mshpoly_t *msh, uint32_t *N_adj);
static void set_nodes_adj_mem(nb_mshpoly_t *msh, uint32_t *N_adj,
			      uint32_t **adj, char *memblock);
static void get_nodes_adj(nb_mshpoly_t *msh, uint32_t *N_adj,
			  uint32_t **adj);
static void set_mask(nb_mshpoly_t *msh, char *mask);

static void update_centroids(nb_mshpoly_t *msh,
			     double (*density)(const double[2],
					       const void *data),
			     const void *density_data);
static void get_centroid(nb_mshpoly_t *msh, uint32_t elem_id, double p[2],
			 double (*density)(const double[2],
					   const void *data),
			 const void *density_data);
static double update_nodes(nb_mshpoly_t *msh, uint32_t *N_adj,
			   uint32_t **adj, char *mask,
			   double (*density)(const double[2],
					     const void *data),
			   const void *density_data);
static void get_new_node_position(nb_mshpoly_t *msh, uint32_t id,
				  uint32_t *N_adj, uint32_t **adj,
				  double new_p[2],
				  double (*density)(const double[2],
					     const void *data),
				  const void *density_data);

uint32_t nb_mshpoly_get_memsize(void)
{
	return sizeof(nb_mshpoly_t);
}

void nb_mshpoly_init(void *mshpoly)
{
	memset(mshpoly, 0, nb_mshpoly_get_memsize());
}

void nb_mshpoly_copy(void *dest, const void *const src)
{
	memcpy(dest, src, nb_mshpoly_get_memsize());
	nb_mshpoly_t *poly = dest;
	const nb_mshpoly_t *const src_poly = src;

	if (poly->N_elems > 0) {
		nb_mshpoly_set_arrays_memory(poly);

		copy_nodes(poly, src_poly);
		copy_edges(poly, src_poly);
		copy_centroids(poly, src_poly);
		copy_N_adj(poly, src_poly);

		nb_mshpoly_set_mem_of_adj_and_ngb(poly);
		copy_adj_and_ngb(poly, src_poly);

		copy_elem_vtx(poly, src_poly);
		copy_N_nod_x_sgm(poly, src_poly);

		nb_mshpoly_set_mem_of_nod_x_sgm(poly);
		copy_nod_x_sgm(poly, src_poly);
	}
}

void nb_mshpoly_set_arrays_memory(nb_mshpoly_t *poly)
{
	uint32_t nod_size = poly->N_nod * 2 * sizeof(*(poly->nod));
	uint32_t edg_size = poly->N_edg * 2 * sizeof(*(poly->edg));
	uint32_t cen_size = poly->N_elems * 2 * sizeof(*(poly->cen));
	uint32_t N_adj_size = poly->N_elems * sizeof(*(poly->N_adj));
	uint32_t adj_size = poly->N_elems * sizeof(*(poly->adj));
	uint32_t ngb_size = poly->N_elems * sizeof(*(poly->ngb));
	uint32_t elem_vtx_size = poly->N_vtx * sizeof(*(poly->vtx));
	uint32_t N_nod_x_sgm_size = poly->N_sgm * sizeof(*(poly->N_nod_x_sgm));
	uint32_t nod_x_sgm_size = poly->N_sgm * sizeof(*(poly->nod_x_sgm));

	uint32_t size = nod_size + edg_size + cen_size + N_adj_size +
		adj_size + ngb_size + elem_vtx_size +
		N_nod_x_sgm_size + nod_x_sgm_size;

	char *memblock = nb_allocate_mem(size);

	poly->nod = (void*) memblock;
	poly->edg = (void*) ((char*)(poly->nod) + nod_size);
	poly->cen = (void*) ((char*)(poly->edg) + edg_size);
	poly->N_adj = (void*) ((char*)(poly->cen) + cen_size);
	poly->adj = (void*) ((char*)(poly->N_adj) + N_adj_size);
	poly->ngb = (void*) ((char*)(poly->adj) + adj_size);
	poly->vtx = (void*) ((char*)(poly->ngb) + ngb_size);
	poly->N_nod_x_sgm = (void*) ((char*)(poly->vtx) +
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

void nb_mshpoly_set_mem_of_adj_and_ngb(nb_mshpoly_t *poly)
{
	uint32_t memsize = get_size_of_adj_and_ngb(poly);
	char *memblock1 = nb_allocate_mem(memsize);
	char *memblock2 = memblock1 + memsize / 2;
	for (uint32_t i = 0; i < poly->N_elems; i++) {
		poly->adj[i] = (void*) memblock1;
		poly->ngb[i] = (void*) memblock2;
		memblock1 += poly->N_adj[i] * sizeof(**(poly->adj));
		memblock2 += poly->N_adj[i] * sizeof(**(poly->ngb));
	}
}

static uint32_t get_size_of_adj_and_ngb(const nb_mshpoly_t *const poly)
{
	uint32_t size = 0;
	for (uint32_t i = 0; i < poly->N_elems; i++)
		size += poly->N_adj[i];
	return 2 * size * sizeof(**(poly->adj));
}

static void copy_adj_and_ngb(nb_mshpoly_t* poly,
			     const nb_mshpoly_t *const src_poly)
{
	for (int i = 0; i < poly->N_elems; i++) {
		memcpy(&(poly->adj[i]), &(src_poly->adj[i]),
		       poly->N_adj[i] * sizeof(**(poly->adj)));
		memcpy(&(poly->ngb[i]), &(src_poly->ngb[i]),
		       poly->N_adj[i] * sizeof(**(poly->ngb)));
	}

}

static void copy_elem_vtx(nb_mshpoly_t* poly,
			  const nb_mshpoly_t *const src_poly)
{
	memcpy(poly->vtx, src_poly->vtx,
	       poly->N_vtx * sizeof(*(poly->vtx)));
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

void nb_mshpoly_set_mem_of_nod_x_sgm(nb_mshpoly_t *poly)
{
	uint32_t memsize = get_size_of_nod_x_sgm(poly);
	char *memblock = nb_allocate_mem(memsize);
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
	nb_mshpoly_t *poly = nb_allocate_mem_poly();
	nb_mshpoly_init(poly);
	return poly;
}

static void* nb_allocate_mem_poly(void)
{
	uint32_t size = nb_mshpoly_get_memsize();
	nb_mshpoly_t *poly = nb_allocate_mem(size);
	return poly;
}

void* nb_mshpoly_clone(const void *const mshpoly_ptr)
{
	nb_mshpoly_t *poly = nb_allocate_mem_poly();
	nb_mshpoly_copy(poly, mshpoly_ptr);
	return poly;
}

void nb_mshpoly_destroy(void *mshpoly_ptr)
{
	nb_mshpoly_finish(mshpoly_ptr);
	nb_free_mem(mshpoly_ptr);
}

void nb_mshpoly_clear(void *mshpoly_ptr)
{
	nb_mshpoly_t *poly = mshpoly_ptr;
	if (NULL != poly->nod) {
		nb_free_mem(poly->adj[0]);
		nb_free_mem(poly->nod_x_sgm[0]);
		nb_free_mem(poly->nod);		
	}
	memset(mshpoly_ptr, 0, nb_mshpoly_get_memsize());
}

uint32_t nb_mshpoly_get_N_invtx(const void *msh)
{
	const nb_mshpoly_t *poly = msh;
	return poly->N_vtx;
}

uint32_t nb_mshpoly_get_N_insgm(const void *msh)
{
	const nb_mshpoly_t *poly = msh;
	return poly->N_sgm;
}

uint32_t nb_mshpoly_get_N_nodes(const void *msh)
{
	const nb_mshpoly_t *poly = msh;
	return poly->N_nod;
}

uint32_t nb_mshpoly_get_N_edges(const void *msh)
{
	const nb_mshpoly_t *poly = msh;
	return poly->N_edg;
}

uint32_t nb_mshpoly_get_N_elems(const void *msh)
{
	const nb_mshpoly_t *poly = msh;
	return poly->N_elems;
}

double nb_mshpoly_node_get_x(const void *msh, uint32_t id)
{
	const nb_mshpoly_t *poly = msh;
	return poly->nod[id * 2];
}

double nb_mshpoly_node_get_y(const void *msh, uint32_t id)
{
	const nb_mshpoly_t *poly = msh;
	return poly->nod[id*2+1];
}

uint32_t nb_mshpoly_edge_get_1n(const void *msh, uint32_t id)
{
	const nb_mshpoly_t *poly = msh;
	return poly->edg[id * 2];
}

uint32_t nb_mshpoly_edge_get_2n(const void *msh, uint32_t id)
{
	const nb_mshpoly_t *poly = msh;
	return poly->edg[id*2+1];
}

void nb_mshpoly_edge_get_midpoint(const void *msh,
				  uint32_t face_id, double w,
				  double midpoint[2])
{
	const nb_mshpoly_t *mshpoly = msh;
	uint32_t n1 = nb_mshpoly_edge_get_1n(msh, face_id);
	uint32_t n2 = nb_mshpoly_edge_get_2n(msh, face_id);
	double *s1 = &(mshpoly->nod[n1 * 2]);
	double *s2 = &(mshpoly->nod[n2 * 2]);

	midpoint[0] = (1 - w) * s1[0] + w * s2[0];
	midpoint[1] = (1 - w) * s1[1] + w * s2[1];
}

double nb_mshpoly_edge_get_normal(const void *msh, uint32_t face_id,
				  double normal[2])
{

	const nb_mshpoly_t *poly = msh;
	uint32_t n1 =  poly->edg[face_id * 2];
	uint32_t n2 =  poly->edg[face_id*2+1];
	double *s1 = &(poly->nod[n1 * 2]);
	double *s2 = &(poly->nod[n2 * 2]);
	double length = nb_utils2D_get_dist(s1, s2);
	normal[0] =  (s2[1] - s1[1]) / length;
	normal[1] = -(s2[0] - s1[0]) / length;
	return length;
}

double nb_mshpoly_elem_get_x(const void *msh, uint32_t id)
{
	const nb_mshpoly_t *poly = msh;
	return poly->cen[id * 2];
}

double nb_mshpoly_elem_get_y(const void *msh, uint32_t id)
{
	const nb_mshpoly_t *poly = msh;
	return poly->cen[id*2+1];
}

double nb_mshpoly_elem_get_area(const void *msh, uint32_t id)
{
	double area = 0.0;

	uint16_t N_adj = nb_mshpoly_elem_get_N_adj(msh, id);
	for (uint16_t i = 0; i < N_adj; i++) {
		uint32_t n1 = nb_mshpoly_elem_get_adj(msh, id, i);
		uint32_t n2 = nb_mshpoly_elem_get_adj(msh, id,
						      (i+1) % N_adj);
		double x1 = nb_mshpoly_node_get_x(msh, n1);
		double y1 = nb_mshpoly_node_get_y(msh, n1);
		double x2 = nb_mshpoly_node_get_x(msh, n2);
		double y2 = nb_mshpoly_node_get_y(msh, n2);
		area += x1 * y2 - x2 * y1;
	}
	return 0.5 * area;
}

double nb_mshpoly_elem_get_radius(const void *msh, uint32_t id)
{
	const nb_mshpoly_t *poly = msh;
	double *x = &(poly->cen[id * 2]);
	uint32_t N_adj = nb_mshpoly_elem_get_N_adj(msh, id);
	double max = 0;
	for (uint16_t i = 0; i < N_adj; i++) {
		uint32_t nid = nb_mshpoly_elem_get_adj(msh, id, i);
		double *ni = &(poly->nod[nid * 2]);
		double dist2 = nb_utils2D_get_dist2(x, ni);
		if (dist2 > max)
			max = dist2;
	}
	return sqrt(max);	
}

double nb_mshpoly_elem_get_apotem(const void *msh, uint32_t id)
{
	const nb_mshpoly_t *poly = msh;
	double *x = &(poly->cen[id * 2]);
	uint32_t N_adj = nb_mshpoly_elem_get_N_adj(msh, id);
	double max = 0;
	for (uint16_t i = 0; i < N_adj; i++) {
		uint32_t id1 = nb_mshpoly_elem_get_adj(msh, id, i);
		uint32_t id2 = nb_mshpoly_elem_get_adj(msh, id, (i+1)%N_adj);
		double *n1 = &(poly->nod[id1 * 2]);
		double *n2 = &(poly->nod[id2 * 2]);
		double closest[2];
		nb_utils2D_get_closest_pnt_to_sgm(n1, n2, x, closest);
		double dist2 = nb_utils2D_get_dist2(x, closest);
		if (dist2 > max)
			max = dist2;
	}
	return sqrt(max);
}

uint32_t nb_mshpoly_elem_find_edge(const void *msh, uint32_t id,
				   uint16_t local_face_id)
{
	const nb_mshpoly_t *poly = msh;
	uint16_t N_adj = poly->N_adj[id];

	uint16_t l1 = local_face_id;
	uint16_t l2 = (local_face_id + 1) % N_adj;

	uint32_t n1 = poly->adj[id][l1];
	uint32_t n2 = poly->adj[id][l2];

	uint32_t N_edges = poly->N_edg;
	uint32_t edge_id = N_edges;
	for (uint32_t i = 0; i < N_edges; i++) {
		uint32_t v1 = poly->edg[i * 2];
		uint32_t v2 = poly->edg[i*2+1];
		if ((n1 == v1 && n2 == v2) || (n1 == v2 && n2 == v1)) {
			edge_id = i;
			break;
		}
	}
	return edge_id;
}

double nb_mshpoly_elem_face_get_length(const void *msh, 
				       uint32_t elem_id,
				       uint16_t face_id)
{
	const nb_mshpoly_t *mshpoly = msh;
	uint16_t N_adj = nb_mshpoly_elem_get_N_adj(msh, elem_id);
	uint32_t n1 = nb_mshpoly_elem_get_adj(msh, elem_id, face_id);
	uint32_t n2 = nb_mshpoly_elem_get_adj(msh, elem_id,
					      (face_id + 1) % N_adj);
	double *s1 = &(mshpoly->nod[n1 * 2]);
	double *s2 = &(mshpoly->nod[n2 * 2]);
	return nb_utils2D_get_dist(s1, s2);
}

double nb_mshpoly_elem_face_get_normal(const void *msh, uint32_t elem_id,
				       uint16_t face_id, double normal[2])
{
	const nb_mshpoly_t *mshpoly = msh;
	uint8_t N_adj = nb_mshpoly_elem_get_N_adj(msh, elem_id);
	uint32_t n1 = nb_mshpoly_elem_get_adj(msh, elem_id, face_id);
	uint32_t n2 = nb_mshpoly_elem_get_adj(msh, elem_id,
					      (face_id + 1) % N_adj);
	double *s1 = &(mshpoly->nod[n1 * 2]);
	double *s2 = &(mshpoly->nod[n2 * 2]);
	double length = nb_utils2D_get_dist(s1, s2);
	normal[0] = (s2[1] - s1[1]) / length;
	normal[1] = -(s2[0] - s1[0]) / length;
	return length;
}

double nb_mshpoly_elem_ngb_get_normal(const void *msh, uint32_t elem_id,
				      uint16_t ngb_id, double normal[2])
{
	const nb_mshpoly_t *mshpoly = msh;
	uint32_t N_elems = nb_mshpoly_get_N_elems(msh);
	uint32_t nid = nb_mshpoly_elem_get_ngb(msh, elem_id, ngb_id);
	memset(normal, 0, 2 * sizeof(double));
	double dist = 0;
	if (nid < N_elems) {
		double *id1 = &(mshpoly->cen[elem_id * 2]);
		double *id2 = &(mshpoly->cen[nid * 2]);

		dist = nb_utils2D_get_dist(id1, id2);
		normal[0] = (id2[0] - id1[0]) / dist;
		normal[1] = (id2[1] - id1[1]) / dist;
	}
	return dist;
}

uint32_t nb_mshpoly_elem_get_N_adj(const void *msh, uint32_t id)
{
	const nb_mshpoly_t *poly = msh;
	return poly->N_adj[id];
}

uint32_t nb_mshpoly_elem_get_adj(const void *msh,
				 uint32_t elem_id, uint8_t adj_id)
{
	const nb_mshpoly_t *poly = msh;
	return poly->adj[elem_id][adj_id];
}

uint32_t nb_mshpoly_elem_get_ngb(const void *msh,
				 uint32_t elem_id, uint8_t ngb_id)
{
	const nb_mshpoly_t *poly = msh;
	return poly->ngb[elem_id][ngb_id];
}

bool nb_mshpoly_elem_has_ngb(const void *msh, uint32_t elem_id,
			     uint16_t ngb_id)
{
	uint32_t N_elems = nb_mshpoly_get_N_elems(msh);
	uint32_t id = nb_mshpoly_elem_get_ngb(msh, elem_id, ngb_id);
	return id < N_elems;
}

bool nb_mshpoly_elem_is_boundary(const void *msh, uint32_t elem_id)
{
	uint16_t i;
	uint16_t N = nb_mshpoly_elem_get_N_adj(msh, elem_id);
	bool out = false;
	for (i = 0; i < N; i++) {
		if (!nb_mshpoly_elem_has_ngb(msh, elem_id, i)) {
			out = true;
			break;
		}
	}
	return out;
}

uint32_t nb_mshpoly_get_invtx(const void *msh, uint32_t id)
{
	const nb_mshpoly_t *poly = msh;
	return poly->vtx[id];
}

uint32_t nb_mshpoly_insgm_get_N_nodes(const void *msh, uint32_t id)
{
	const nb_mshpoly_t *poly = msh;
	return poly->N_nod_x_sgm[id];
}

uint32_t nb_mshpoly_insgm_get_node(const void *msh, uint32_t sgm_id,
				     uint32_t node_id)
{
	const nb_mshpoly_t *poly = msh;
	return poly->nod_x_sgm[sgm_id][node_id];
}

void nb_mshpoly_get_enveloping_box(const void *msh, double box[4])
{
	const nb_mshpoly_t *poly = msh;
	nb_utils2D_get_enveloping_box(poly->N_nod, poly->nod,
				       2 * sizeof(*(poly->nod)),
				       nb_utils2D_get_x_from_darray,
				       nb_utils2D_get_y_from_darray,
				       box);
}

bool nb_mshpoly_is_vtx_inside(const void *msh, double x, double y)
{
	return false;/* PENDING */
}

void nb_mshpoly_build_model(const void *msh, nb_model_t *model)
{
	/* PENDING */
}

void nb_mshpoly_build_model_disabled_elems(const void *msh,
					   const bool *elems_enabled,
					   nb_model_t *model,
					   uint32_t *N_input_vtx,
					   uint32_t **input_vtx)
{
	/* PENDING */
}

void nb_mshpoly_load_from_tessellator2D(void *mshpoly, nb_tessellator2D_t *mesh)
{
	if (nb_tessellator2D_get_N_trg(mesh) > 0) {
		split_exterior_trg(mesh);

		mesh_enumerate_vtx(mesh);
		mesh_enumerate_trg(mesh);

		vinfo_t vinfo;
		init_voronoi_info(&vinfo, mesh);

		uint32_t Nt = nb_tessellator2D_get_N_trg(mesh);
		uint32_t *trg_cc_map = nb_allocate_mem(Nt * sizeof(uint32_t));
		init_trg_cc_map(trg_cc_map, Nt);

		vgraph_t vgraph;
		init_voronoi_graph(&vgraph, &vinfo, trg_cc_map, mesh);

		create_mapping(&vinfo, &vgraph, mesh, trg_cc_map);
		nb_free_mem(trg_cc_map);

		set_voronoi(mshpoly, &vgraph, &vinfo, mesh);

		finish_voronoi_info(&vinfo);
		finish_voronoi_graph(&vgraph);
	}
}

static void init_voronoi_info(vinfo_t *vinfo,
			      const nb_tessellator2D_t *const mesh)
{
	memset(vinfo, 0,  sizeof(*vinfo));
	uint32_t Nv = nb_tessellator2D_get_N_vtx(mesh);
	uint32_t Nt = nb_tessellator2D_get_N_trg(mesh);
	uint32_t size1 = Nv * sizeof(uint32_t);
	uint32_t size2 = Nt * sizeof(uint32_t);
	char *memblock = nb_allocate_mem(2 * size1 + size2);
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
			       const nb_tessellator2D_t *const mesh)
{
	vgraph->N = nb_tessellator2D_get_N_vtx(mesh);
	
	uint32_t N_edg = nb_tessellator2D_get_N_edg(mesh);

	uint32_t size1 = vgraph->N * sizeof(*(vgraph->N_adj));
	uint32_t size2 = vgraph->N * sizeof(*(vgraph->adj));
	uint32_t memsize = size1 + size2 +
		2 * N_edg * sizeof(**(vgraph->adj));
	char *memblock = nb_allocate_mem(memsize);

	vgraph->N_adj = (void*) (memblock);
	vgraph->adj = (void*) (memblock +  size1);

	count_vgraph_adj(vgraph, mesh);
	set_vgraph_adj_mem(vgraph, memblock + size1 + size2);
	set_vgraph_adj(vgraph, vinfo, trg_cc_map, mesh);
}

static void create_mapping(vinfo_t *vinfo,
			   const vgraph_t *const vgraph,
			   const nb_tessellator2D_t *const mesh,
			   const uint32_t *trg_cc_map)
{
	uint32_t ielem = 0;
	uint32_t inode = 0;

	nb_bins2D_iter_t* biter = nb_allocate_on_stack(nb_bins2D_iter_get_memsize());
	nb_bins2D_iter_init(biter);
	nb_bins2D_iter_set_bins(biter, mesh->ug_vtx);
	while (nb_bins2D_iter_has_more(biter)) {
		const msh_vtx_t* vtx = nb_bins2D_iter_get_next(biter);
		uint32_t id = mvtx_get_id(vtx);
		if (mvtx_is_type_location(vtx, INTERIOR)) {
			vinfo->N_vtx_in += 1;
			vinfo->vtx_map[id] = ielem;
			ielem += 1;
		} else {
			vinfo->N_vtx_out += 1;
			vinfo->vtx_map[id] = inode;
			inode += 1;
		}
	}
	nb_bins2D_iter_finish(biter);

	uint16_t iter_size = nb_iterator_get_memsize();
	nb_iterator_t *iter = nb_allocate_on_stack(iter_size);
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, mesh->ht_trg);
	while (nb_iterator_has_more(iter)) {
		const msh_trg_t *trg = nb_iterator_get_next(iter);
		uint32_t id = trg->id;
		if (trg_is_interior(trg, vgraph)) {
			if (id == trg_cc_map[id]) {
				vinfo->trg_map[id] = inode;
				inode += 1;
			}
		}
	}
	nb_iterator_restart(iter);
	while (nb_iterator_has_more(iter)) {
		const msh_trg_t *trg = nb_iterator_get_next(iter);
		uint32_t id = trg->id;
		if (trg_is_interior(trg, vgraph)) {
			if (id != trg_cc_map[id]) {
				uint32_t cc_id = trg_cc_map[id];
				vinfo->trg_map[id] = vinfo->trg_map[cc_id];
			}
		}
	}
	nb_iterator_finish(iter);
	
}

static bool trg_is_interior(const msh_trg_t *const trg,
			    const vgraph_t *const vgraph)
{
	return mvtx_is_type_location(trg->v1, INTERIOR) &&
		mvtx_is_type_location(trg->v2, INTERIOR) &&
		mvtx_is_type_location(trg->v3, INTERIOR);

}

static void count_vgraph_adj(vgraph_t *vgraph,
			     const nb_tessellator2D_t *const mesh)
{
	memset(vgraph->N_adj, 0,  vgraph->N * sizeof(*(vgraph->N_adj)));

	uint16_t iter_size = nb_iterator_get_memsize();
	nb_iterator_t *iter = nb_allocate_on_stack(iter_size);
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, mesh->ht_edge);
	while (nb_iterator_has_more(iter)) {
		const msh_edge_t *edge = nb_iterator_get_next(iter);
		uint32_t id1 = mvtx_get_id(edge->v1);
		uint32_t id2 = mvtx_get_id(edge->v2);

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

static void set_vgraph_adj(vgraph_t *vgraph, vinfo_t *vinfo,
			   uint32_t *trg_cc_map,
			   const nb_tessellator2D_t *const mesh)
{
	memset(vgraph->N_adj, 0,  vgraph->N * sizeof(*(vgraph->N_adj)));

	uint16_t iter_size = nb_iterator_get_memsize();
	nb_iterator_t *iter = nb_allocate_on_stack(iter_size);
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, mesh->ht_edge);
	while (nb_iterator_has_more(iter)) {
		const msh_edge_t *edge = nb_iterator_get_next(iter);

		bool cc = adj_is_cocircular(edge);
		update_cc_map(edge, cc, trg_cc_map,
			      nb_tessellator2D_get_N_trg(mesh));
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
		msh_vtx_t *v4 = mtrg_get_opposite_vertex(edg->t2, edg);
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
		uint32_t id1 = edge->t1->id;
		uint32_t id2 = edge->t2->id;

		uint32_t new_cc = trg_cc_map[id1];
		uint32_t old_cc = trg_cc_map[id2];
		for (uint32_t i = 0; i < N_trg; i++) {
			if (old_cc == trg_cc_map[i])
				trg_cc_map[i] = new_cc;
		}
	}
}

static void insert_edg_as_adj(vgraph_t *vgraph, const msh_edge_t *edge,
			      bool is_cc)
{
	bool semi_onsegment =
		mvtx_is_type_location(edge->v1, ONSEGMENT) ||
		mvtx_is_type_location(edge->v2, ONSEGMENT);
	bool not_interior_cocircular = !is_cc || semi_onsegment;
	bool semi_interior =
		mvtx_is_type_location(edge->v1, INTERIOR) ||
		mvtx_is_type_location(edge->v2, INTERIOR);
	if (not_interior_cocircular && semi_interior) {
		insert_adj_sorted_by_angle(vgraph, edge, edge->v1, edge->v2);
		insert_adj_sorted_by_angle(vgraph, edge, edge->v2, edge->v1);
	}
}

static void insert_adj_sorted_by_angle(vgraph_t *vgraph,
				       const msh_edge_t *edge,
				       const msh_vtx_t *v1,
				       const msh_vtx_t *v2)
{
	uint32_t id_global = mvtx_get_id(v1);
	uint32_t id_local = vgraph->N_adj[id_global];
	double angle_id;
	if (0 < id_local) {
		double x = v2->x[0] - v1->x[0];
		double y = v2->x[1] - v1->x[1];
		angle_id = atan2(y, x);
	}

	for (uint16_t j = 0; j < id_local; j++) {
		msh_edge_t *jedge = vgraph->adj[id_global][j];
		v2 = medge_get_partner_vtx(jedge, v1);
		double x = v2->x[0] - v1->x[0];
		double y = v2->x[1] - v1->x[1];
		double angle_j = atan2(y, x);
		if (angle_j > angle_id) {
			msh_edge_t *aux = vgraph->adj[id_global][j];
			vgraph->adj[id_global][j] = (msh_edge_t*) edge;
			edge = aux;
			angle_id = angle_j;
		}
	}
	vgraph->adj[id_global][id_local] = (msh_edge_t*) edge;
	vgraph->N_adj[id_global] += 1;
}

static void counting_edg_in_vinfo(vinfo_t *vinfo, const vgraph_t *vgraph,
				  const msh_edge_t *edge, bool is_cc)
{
	bool extremes_onsegment =
		mvtx_is_type_location(edge->v1, ONSEGMENT) &&
		mvtx_is_type_location(edge->v2, ONSEGMENT);
	bool edge_is_interior =
		mvtx_is_type_location(edge->v1, INTERIOR) &&
		mvtx_is_type_location(edge->v2, INTERIOR);
	if (extremes_onsegment) {
		vinfo->N_edg_out += 1;
	} else if (edge_is_interior) {
		vinfo->N_edg_in += 1;

		msh_vtx_t *opp1 = mtrg_get_opposite_vertex(edge->t1, edge);
		msh_vtx_t *opp2 = mtrg_get_opposite_vertex(edge->t2, edge);
		bool t1_is_interior = mvtx_is_type_location(opp1, INTERIOR);
		bool t2_is_interior = mvtx_is_type_location(opp2, INTERIOR);
		if (is_cc && t1_is_interior && t2_is_interior)
			vinfo->N_cc_in += 1;

		if (t1_is_interior)
			vinfo->N_trg_in += 1;

		if (t2_is_interior)
			vinfo->N_trg_in += 1;
	}
}

static void finish_voronoi_info(vinfo_t *vinfo)
{
	nb_free_mem(vinfo->vtx_map);
}

static void finish_voronoi_graph(vgraph_t *vgraph)
{
	nb_free_mem(vgraph->N_adj);
}

static void set_voronoi(nb_mshpoly_t *poly,
			const vgraph_t *const vgraph,
			const vinfo_t *const vinfo,
			const nb_tessellator2D_t *const mesh)
{
	set_quantities(poly, vinfo, mesh);
	
	nb_mshpoly_set_arrays_memory(poly);

	set_nodes_and_centroids(poly, vgraph, vinfo, mesh);
	set_edges(poly, vgraph, vinfo, mesh);
	set_N_adj(poly, vgraph, vinfo, mesh);
	
	nb_mshpoly_set_mem_of_adj_and_ngb(poly);
	set_adj_and_ngb(poly, vgraph, vinfo, mesh);

	set_elem_vtx(poly, vgraph, vinfo, mesh);
	set_N_nod_x_sgm(poly, mesh);

	nb_mshpoly_set_mem_of_nod_x_sgm(poly);
	set_nod_x_sgm(poly, vinfo, mesh);
}


static void set_quantities(nb_mshpoly_t *poly,
			   const vinfo_t *const vinfo,
			   const nb_tessellator2D_t *const mesh)
{
	poly->N_nod = vinfo->N_trg_in + vinfo->N_vtx_out - vinfo->N_cc_in;
	poly->N_edg = vinfo->N_edg_in + vinfo->N_edg_out - vinfo->N_cc_in;
	poly->N_elems = vinfo->N_vtx_in;

	poly->N_vtx = mesh->N_input_vtx;
	poly->N_sgm = mesh->N_input_sgm;
}

static void set_nodes_and_centroids(nb_mshpoly_t *poly,
				    const vgraph_t *const vgraph,
				    const vinfo_t *const vinfo,
				    const nb_tessellator2D_t *const mesh)
{
	nb_bins2D_iter_t* biter = nb_allocate_on_stack(nb_bins2D_iter_get_memsize());
	nb_bins2D_iter_init(biter);
	nb_bins2D_iter_set_bins(biter, mesh->ug_vtx);
	while (nb_bins2D_iter_has_more(biter)) {
		msh_vtx_t* vtx = (void*) nb_bins2D_iter_get_next(biter);
		uint32_t id = mvtx_get_id(vtx);
		if (mvtx_is_type_location(vtx, INTERIOR)) {
			uint32_t ielem = vinfo->vtx_map[id];
			scale_vtx(&(poly->cen[ielem * 2]), vtx->x, mesh);
		} else {
			/* TEMPORAL: Enhance estimation at boundaries */
			uint32_t inode = vinfo->vtx_map[id];
			scale_vtx(&(poly->nod[inode * 2]), vtx->x, mesh);
		}
	}
	nb_bins2D_iter_finish(biter);

	uint16_t iter_size = nb_iterator_get_memsize();
	nb_iterator_t *iter = nb_allocate_on_stack(iter_size);
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, mesh->ht_trg);
	while (nb_iterator_has_more(iter)) {
		const msh_trg_t *trg = nb_iterator_get_next(iter);
		uint32_t id = trg->id;
		if (trg_is_interior(trg, vgraph)) {
			double circumcenter[2];
			nb_utils2D_get_circumcenter(trg->v1->x,
						     trg->v2->x,
						     trg->v3->x,
						     circumcenter);
			uint32_t inode = vinfo->trg_map[id];
			scale_vtx(&(poly->nod[inode * 2]),
				  circumcenter, mesh);;
		}
	}
	nb_iterator_finish(iter);
}

static void scale_vtx(double out[2], double in[2],
		      const nb_tessellator2D_t *mesh)
{
	out[0] = in[0] / mesh->scale + mesh->xdisp;
	out[1] = in[1] / mesh->scale + mesh->ydisp;
}

static void set_edges(nb_mshpoly_t *poly,
		      const vgraph_t *const vgraph,
		      const vinfo_t *const vinfo,
		      const nb_tessellator2D_t *const mesh)
{
	uint32_t iedge = 0;

	uint16_t iter_size = nb_iterator_get_memsize();
	nb_iterator_t *iter = nb_allocate_on_stack(iter_size);
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, mesh->ht_edge);
	while (nb_iterator_has_more(iter)) {
		const msh_edge_t *edge = nb_iterator_get_next(iter);
		uint32_t v1 = mvtx_get_id(edge->v1);
		uint32_t v2 = mvtx_get_id(edge->v2);
		bool extremes_onsegment =
			mvtx_is_type_location(edge->v1, ONSEGMENT) &&
			mvtx_is_type_location(edge->v2, ONSEGMENT);
		bool edge_is_interior =
			mvtx_is_type_location(edge->v1, INTERIOR) &&
			mvtx_is_type_location(edge->v2, INTERIOR);
		if (extremes_onsegment) {
			/* Process not interior edges */
			poly->edg[iedge * 2] = vinfo->vtx_map[v1];
			poly->edg[iedge*2+1] = vinfo->vtx_map[v2];
			iedge += 1;
		} else if (edge_is_interior) {
			process_interior_edge(poly, vgraph, vinfo,
					      edge, iedge);
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
	uint32_t opp1 = mvtx_get_id(opp_t1);
	uint32_t opp2 = mvtx_get_id(opp_t2);
	bool t1_is_interior = mvtx_is_type_location(opp_t1, INTERIOR);
	bool t2_is_interior = mvtx_is_type_location(opp_t2, INTERIOR);

	uint32_t t1 = edg->t1->id;
	uint32_t t2 = edg->t2->id;
	bool is_not_cc = vinfo->trg_map[t1] != vinfo->trg_map[t2];
	if (is_not_cc && t1_is_interior && t2_is_interior) {
		/* Interior triangles (No cocircular) */
		poly->edg[iedge * 2] = vinfo->trg_map[t1];
		poly->edg[iedge*2+1] = vinfo->trg_map[t2];		
	} else {
		/* At least a triangle is not interior */
		uint32_t node1;
		if (t1_is_interior)
			node1 = vinfo->trg_map[t1];
		else
			node1 = vinfo->vtx_map[opp1];
		uint32_t node2;
		if (t2_is_interior)
			node2 = vinfo->trg_map[t2];
		else
			node2 = vinfo->vtx_map[opp2];

		poly->edg[iedge * 2] = node1;
		poly->edg[iedge*2+1] = node2;			
	}
}

static void set_N_adj(nb_mshpoly_t *poly,
		      const vgraph_t *const vgraph,
		      const vinfo_t *const vinfo,
		      const nb_tessellator2D_t *const mesh)
{
	nb_bins2D_iter_t* biter = nb_allocate_on_stack(nb_bins2D_iter_get_memsize());
	nb_bins2D_iter_init(biter);
	nb_bins2D_iter_set_bins(biter, mesh->ug_vtx);
	while (nb_bins2D_iter_has_more(biter)) {
		const msh_vtx_t* vtx = nb_bins2D_iter_get_next(biter);
		if (mvtx_is_type_location(vtx, INTERIOR)) {
			uint32_t id = mvtx_get_id(vtx);
			uint32_t elem_id = vinfo->vtx_map[id];
			poly->N_adj[elem_id] = vgraph->N_adj[id];
		}
	}
	nb_bins2D_iter_finish(biter);
}

static void set_adj_and_ngb(nb_mshpoly_t *poly,
			    const vgraph_t *const vgraph,
			    const vinfo_t *const vinfo,
			    const nb_tessellator2D_t *const mesh)
{
	nb_bins2D_iter_t* biter = nb_allocate_on_stack(nb_bins2D_iter_get_memsize());
	nb_bins2D_iter_init(biter);
	nb_bins2D_iter_set_bins(biter, mesh->ug_vtx);
	while (nb_bins2D_iter_has_more(biter)) {
		const msh_vtx_t* vtx = nb_bins2D_iter_get_next(biter);
		if (mvtx_is_type_location(vtx, INTERIOR)) {
			uint32_t id = mvtx_get_id(vtx);
			uint16_t id_adj = 0;
			for (uint16_t j = 0; j < vgraph->N_adj[id]; j++) {
				id_adj = add_adj_and_ngb(poly, vgraph,
							 vinfo, id, j,
							 id_adj);
			}
			uint32_t elem_id = vinfo->vtx_map[id];
			poly->N_adj[elem_id] = id_adj;
		}
	}
	nb_bins2D_iter_finish(biter);
}

static uint16_t add_adj_and_ngb(nb_mshpoly_t *poly,
				const vgraph_t *const vgraph,
				const vinfo_t *const vinfo,
				uint32_t i, uint16_t j, uint16_t id_adj)
{
	uint32_t elem_id = vinfo->vtx_map[i];
	uint16_t id_prev = (j + vgraph->N_adj[i] - 1) % vgraph->N_adj[i];
	msh_vtx_t *v1 = get_partner(vgraph, vinfo, i, id_prev);
	msh_vtx_t *v2 = get_partner(vgraph, vinfo, i, j);
	uint32_t id1 = mvtx_get_id(v1);
	uint32_t id2 = mvtx_get_id(v2);
	bool v1_is_interior = mvtx_is_type_location(v1, INTERIOR);
	bool v2_is_interior = mvtx_is_type_location(v2, INTERIOR);
	if (v2_is_interior) {
		if (v1_is_interior) {
			/* Interior trg case */
			msh_trg_t *trg = get_prev_trg(vgraph, vinfo, i, j);
			uint32_t trg_id = trg->id;
			poly->adj[elem_id][id_adj] = vinfo->trg_map[trg_id];
		} else {
			/* First node in the boundary */
			poly->adj[elem_id][id_adj] = vinfo->vtx_map[id1];
		}
		poly->ngb[elem_id][id_adj] = vinfo->vtx_map[id2];
		id_adj += 1;
	} else {
		if (!v1_is_interior) {
			/* Both nodes are in the boundary */
			poly->adj[elem_id][id_adj] = vinfo->vtx_map[id1];
			poly->ngb[elem_id][id_adj] = poly->N_elems;
			id_adj += 1;
		}
	}
	return id_adj;
}

static msh_vtx_t *get_partner(const vgraph_t *const vgraph,
			      const vinfo_t *const vinfo,
			      uint32_t i, uint16_t j)
{
	msh_vtx_t *out = NULL;
	if (j < vgraph->N_adj[i]) {
		msh_edge_t *edge = vgraph->adj[i][j];
		if (NULL != edge) {
			if (i == mvtx_get_id(edge->v1))
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
	msh_trg_t *trg;
	msh_edge_t *edge = vgraph->adj[i][j];
	if (i == mvtx_get_id(edge->v1))
		trg = edge->t2;
	else if (i == mvtx_get_id(edge->v2))
		trg = edge->t1;
	else
		trg = NULL;
	return trg;
}

static void set_elem_vtx(nb_mshpoly_t *poly,
			 const vgraph_t *const vgraph,
			 const vinfo_t *const vinfo,
			 const nb_tessellator2D_t *const mesh)
{
	for (uint32_t i = 0; i < poly->N_vtx; i++) {
		msh_vtx_t *vtx = mesh->input_vtx[i];
		if (mvtx_is_type_location(vtx, ONSEGMENT)) {
			int id = mvtx_get_id(vtx);
			poly->vtx[i] = vinfo->vtx_map[id];
		} else {
			poly->vtx[i] = poly->N_elems;
		}
	}
}

static void set_N_nod_x_sgm(nb_mshpoly_t *poly,
			    const nb_tessellator2D_t *const mesh)
{
	for (uint32_t i = 0; i < poly->N_sgm; i++) {
		msh_edge_t* sgm = mesh->input_sgm[i];
		uint32_t counter = 0;
		while (NULL != sgm) {
			counter++;
			sgm = medge_subsgm_next(sgm);
		}
		poly->N_nod_x_sgm[i] = counter + 1;
	}

}

static void set_nod_x_sgm(nb_mshpoly_t *poly,
			  const vinfo_t *const vinfo,
			  const nb_tessellator2D_t *const mesh)
{
	for (uint32_t i = 0; i < poly->N_sgm; i++)
		set_sgm_nodes(poly, vinfo, mesh, i);
}

static void set_sgm_nodes(nb_mshpoly_t *poly,
			  const vinfo_t *const vinfo,
			  const nb_tessellator2D_t *const mesh,
			  uint32_t sgm_id)
{
	msh_edge_t *sgm_prev = mesh->input_sgm[sgm_id];
	if (NULL != sgm_prev) {
		msh_edge_t *sgm = medge_subsgm_next(sgm_prev);
		if (NULL == sgm) {
			uint32_t id1 = mvtx_get_id(sgm_prev->v1);
			uint32_t id2 = mvtx_get_id(sgm_prev->v2);
			poly->nod_x_sgm[sgm_id][0] = vinfo->vtx_map[id1];
			poly->nod_x_sgm[sgm_id][1] = vinfo->vtx_map[id2];
		} else {
			assemble_sgm_wire(poly, vinfo,
					  sgm_id, sgm_prev, sgm);
		}
	}
}

static void assemble_sgm_wire(nb_mshpoly_t *poly,
			      const vinfo_t *const vinfo,
			      uint32_t sgm_id,
			      msh_edge_t *sgm_prev, msh_edge_t *sgm)
{
	uint32_t idx = 0;
	uint32_t id_chain;
	uint32_t id1 = mvtx_get_id(sgm_prev->v1);
	uint32_t id2 = mvtx_get_id(sgm_prev->v2);
	uint32_t id1n = mvtx_get_id(sgm->v1);
	uint32_t id2n = mvtx_get_id(sgm->v2);
	if (id2 == id1n || id2 == id2n) {
		poly->nod_x_sgm[sgm_id][idx++] =  vinfo->vtx_map[id1];
		poly->nod_x_sgm[sgm_id][idx++] =  vinfo->vtx_map[id2];
		id_chain = id2;
	} else {
		poly->nod_x_sgm[sgm_id][idx++] =  vinfo->vtx_map[id2];
		poly->nod_x_sgm[sgm_id][idx++] =  vinfo->vtx_map[id1];
		id_chain = id1;
	}
	while (NULL != sgm) {
		sgm_prev = sgm;
		uint32_t id1 = mvtx_get_id(sgm_prev->v1);
		uint32_t id2 = mvtx_get_id(sgm_prev->v2);
		if (id1 == id_chain) {
			poly->nod_x_sgm[sgm_id][idx++] = vinfo->vtx_map[id2];
			id_chain = id2;
		} else {
			poly->nod_x_sgm[sgm_id][idx++] = vinfo->vtx_map[id1];
			id_chain = id1;
		}
		sgm = medge_subsgm_next(sgm);
	}
}

static void split_exterior_trg(nb_tessellator2D_t *mesh)
{
	nb_container_t *exterior_trg =
		nb_allocate_on_stack(nb_container_get_memsize(NB_QUEUE));
	nb_container_init(exterior_trg, NB_QUEUE);
	
	initialize_exterior_trg(mesh, exterior_trg);
	delete_exterior_trg(mesh, exterior_trg);

	nb_container_finish(exterior_trg);
}

static void initialize_exterior_trg(const nb_tessellator2D_t *mesh,
				    nb_container_t *exterior_trg)
{
	uint32_t size = nb_iterator_get_memsize();
	nb_iterator_t *iter = nb_allocate_on_stack(size);
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, mesh->ht_trg);
	while (nb_iterator_has_more(iter)) {
		const msh_trg_t *trg = nb_iterator_get_next(iter);
		if (is_exterior(trg))
			nb_container_insert(exterior_trg, trg);
	}
	nb_iterator_finish(iter);
}

static bool is_exterior(const msh_trg_t *trg)
{
	/* bool out = mtrg_has_an_input_sgm(trg);
	   if (!out) TEMPORAL commented because 
	   error on boundary */
	bool out = have_all_nodes_in_sgm(trg);
	return out;
}

static bool have_all_nodes_in_sgm(const msh_trg_t *trg)
{
	bool out = false;
	if (mvtx_is_type_location(trg->v1, ONSEGMENT)) {
		if (mvtx_is_type_location(trg->v2, ONSEGMENT)) {
			if (mvtx_is_type_location(trg->v3, ONSEGMENT))
				out = true;
		}
	}
	return out;
}


static void delete_exterior_trg(nb_tessellator2D_t *mesh,
				nb_container_t *exterior_trg)
{
	while (nb_container_is_not_empty(exterior_trg)) {
		msh_trg_t *trg = nb_container_delete_first(exterior_trg);
		if (nb_container_exist(mesh->ht_trg, trg)) {
			msh_vtx_t *cen = mvtx_create(mesh);
			nb_utils2D_trg_get_centroid(trg->v1->x,
						     trg->v2->x,
						     trg->v3->x,
						     cen->x);
			nb_ruppert_insert_verified_vtx(mesh, trg,
						       cen);
		}
	}

}
void nb_mshpoly_set_nodal_permutation(void *msh, const uint32_t *perm)
{
	set_nodal_perm_to_nodes(msh, perm);
	set_nodal_perm_to_edges(msh, perm);
	set_nodal_perm_to_elems(msh, perm);
	set_nodal_perm_to_invtx(msh, perm);
	set_nodal_perm_to_insgm(msh, perm);
}

static void set_nodal_perm_to_nodes(nb_mshpoly_t *poly, const uint32_t *perm)
{
	uint32_t N = poly->N_nod;
	
	uint32_t memsize = N * 2 * sizeof(*(poly->nod));
	double *nodes = nb_soft_allocate_mem(memsize);

	memcpy(nodes, poly->nod, memsize);

	for (uint32_t i = 0; i < N; i++) {
		uint32_t id = perm[i];
		memcpy(&(poly->nod[id*2]), &(nodes[i*2]), 2 * sizeof(*nodes));
	}

	nb_soft_free_mem(memsize, nodes);
}

static void set_nodal_perm_to_edges(nb_mshpoly_t *poly, const uint32_t *perm)
{
	uint32_t N = poly->N_edg;
	for (uint32_t i = 0; i < 2 * N; i++) {
		uint32_t id = poly->edg[i];
		poly->edg[i] = perm[id];
	}
}

static void set_nodal_perm_to_elems(nb_mshpoly_t *poly, const uint32_t *perm)
{
	uint32_t N = poly->N_elems;
	for (uint32_t i = 0; i < N; i++) {
		for (uint16_t j = 0; j < poly->N_adj[i]; j++) {
			uint32_t id = poly->adj[i][j];
			poly->adj[i][j] = perm[id];
		}
	}
}

static void set_nodal_perm_to_invtx(nb_mshpoly_t *poly, const uint32_t *perm)
{
	uint32_t N = poly->N_vtx;
	for (uint32_t i = 0; i < N; i++) {
		uint32_t id = poly->vtx[i];
		poly->vtx[i] = perm[id];
	}
}

static void set_nodal_perm_to_insgm(nb_mshpoly_t *poly, const uint32_t *perm)
{
	uint32_t N = poly->N_sgm;
	for (uint32_t i = 0; i < N; i++) {
		for (uint16_t j = 0; j < poly->N_nod_x_sgm[i]; j++) {
			uint32_t id = poly->nod_x_sgm[i][j];
			poly->nod_x_sgm[i][j] = perm[id];
		}
	}
}

void nb_mshpoly_centroid_iteration(void *mshpoly, uint32_t max_iter,
				   /* density can be NULL */
				   double (*density)(const double[2],
						     const void *data),
				   const void *density_data)
{
	uint32_t N_nodes = nb_mshpoly_get_N_nodes(mshpoly);
	char *memblock1 = nb_allocate_mem(N_nodes * (sizeof(uint32_t) +
						     sizeof(uint32_t*) + 1));
	uint32_t *N_adj = (void*) memblock1;
	uint32_t **adj = (void*) (memblock1 + N_nodes * sizeof(uint32_t));
	                     /* Mask is for fixed nodes (at the boundary) */
	char *mask = (void*) (memblock1 + N_nodes * (sizeof(uint32_t) +
						     sizeof(uint32_t*)));
	uint32_t N_total_adj = get_nodes_N_adj(mshpoly, N_adj);
	
	char *memblock2 = nb_allocate_mem(N_total_adj * sizeof(uint32_t));
	set_nodes_adj_mem(mshpoly, N_adj, adj, memblock2);
	get_nodes_adj(mshpoly, N_adj, adj);
	set_mask(mshpoly, mask);

	double max_disp2 = 1;
	uint32_t i = 0;
	while (max_disp2 > 1e-25 && i < max_iter) {
		update_centroids(mshpoly, density, density_data);
		max_disp2 = update_nodes(mshpoly, N_adj, adj, mask,
					 density, density_data);
		i++;
	}

	nb_free_mem(memblock1);
	nb_free_mem(memblock2);
}

static uint32_t get_nodes_N_adj(nb_mshpoly_t *msh, uint32_t *N_adj)
{
	uint32_t N_nodes = nb_mshpoly_get_N_nodes(msh);
	memset(N_adj, 0, N_nodes * sizeof(*N_adj));

	uint32_t N_elems = nb_mshpoly_get_N_elems(msh);
	uint32_t N_total = 0;
	for (uint32_t i = 0; i < N_elems; i++) {
		int N = nb_mshpoly_elem_get_N_adj(msh, i);
		for (int j = 0; j < N; j++) {
			uint32_t id = nb_mshpoly_elem_get_adj(msh, i, j);
			N_adj[id] += 1;
		}
		N_total += N;
	}
	return N_total;
}

static void set_nodes_adj_mem(nb_mshpoly_t *msh, uint32_t *N_adj,
			      uint32_t **adj, char *memblock)
{
	uint32_t N = nb_mshpoly_get_N_nodes(msh);

	for (uint32_t i = 0; i < N; i++) {
		adj[i] = (void*) memblock;
		memblock += N_adj[i] * sizeof(uint32_t);
	}
}

static void get_nodes_adj(nb_mshpoly_t *msh, uint32_t *N_adj,
			  uint32_t **adj)
{
	uint32_t N_nodes = nb_mshpoly_get_N_nodes(msh);
	memset(N_adj, 0, N_nodes * sizeof(*N_adj));

	uint32_t N_elems = nb_mshpoly_get_N_elems(msh);
	for (uint32_t i = 0; i < N_elems; i++) {
		int N = nb_mshpoly_elem_get_N_adj(msh, i);
		for (int j = 0; j < N; j++) {
			uint32_t id = nb_mshpoly_elem_get_adj(msh, i, j);
			uint32_t local_id = N_adj[id];
			adj[id][local_id] = i;
			N_adj[id] = local_id + 1;
		}
	}
}

static void set_mask(nb_mshpoly_t *msh, char *mask)
/*
 * 0: Free vertices (steiner points).
 * 1: Vertices constrained to input segments.
 * 2: Input vertices fixed.
 */
{
	uint32_t N = nb_mshpoly_get_N_nodes(msh);
	memset(mask, 0, N);

	uint32_t N_sgm = nb_mshpoly_get_N_insgm(msh);
	for (uint32_t i = 0; i < N_sgm; i++) {
		uint32_t N_sgm_nodes = nb_mshpoly_insgm_get_N_nodes(msh, i);
		for (uint32_t j = 0; j < N_sgm_nodes; j++) {
			uint32_t id = nb_mshpoly_insgm_get_node(msh, i, j);
			if (id < N)
				mask[id] = 1;
		}
	}

	uint32_t N_vtx = nb_mshpoly_get_N_invtx(msh);
	for (uint32_t i = 0; i < N_vtx; i++) {
		uint32_t id = nb_mshpoly_get_invtx(msh, i);
		if (id < N)
			mask[id] = 2;
	}
}

static void update_centroids(nb_mshpoly_t *msh,
			     double (*density)(const double[2],
					       const void *data),
			     const void *density_data)
{
	
	uint32_t N_elems = nb_mshpoly_get_N_elems(msh);
	for (uint32_t i = 0; i < N_elems; i++) {
		double p[2];
		get_centroid(msh, i, p, density, density_data);
		msh->cen[i * 2] = p[0];
		msh->cen[i*2+1] = p[1];
	}
}

static void get_centroid(nb_mshpoly_t *msh, uint32_t elem_id, double p[2],
			 double (*density)(const double[2],
					   const void *data),
			 const void *density_data)
{
	double t1[2];
	double t2[2];
	double t3[2];
	uint32_t id1 = nb_mshpoly_elem_get_adj(msh, elem_id, 0);
	t1[0] = nb_mshpoly_node_get_x(msh, id1);
	t1[1] = nb_mshpoly_node_get_y(msh, id1);
	double den1 = 1.0;
	if (NULL != density)
		den1 = density(t1, density_data);

	double total_area = 0.0;
	p[0] = 0.0;
	p[1] = 0.0;

	int N = nb_mshpoly_elem_get_N_adj(msh, elem_id);
	for (int j = 2; j < N; j++) {
		uint32_t id2 = nb_mshpoly_elem_get_adj(msh, elem_id, j-1);
		t2[0] = nb_mshpoly_node_get_x(msh, id2);
		t2[1] = nb_mshpoly_node_get_y(msh, id2);
		double den2 = 1.0;
		if (NULL != density)
			den2 = density(t2, density_data);

		uint32_t id3 = nb_mshpoly_elem_get_adj(msh, elem_id, j);
		t3[0] = nb_mshpoly_node_get_x(msh, id3);
		t3[1] = nb_mshpoly_node_get_y(msh, id3);
		double den3 = 1.0;
		if (NULL != density)
			den3 = density(t3, density_data);

		double ct[2];
		ct[0] = (den1 * t1[0] + den2 * t2[0] + den3 * t3[0]) / 3;
		ct[1] = (den1 * t1[1] + den2 * t2[1] + den3 * t3[1]) / 3;

		double area = nb_utils2D_get_trg_area(t1, t2, t3);
		total_area += area;

		p[0] += ct[0] * area;
		p[1] += ct[1] * area;
	}
	p[0] /= total_area;
	p[1] /= total_area;
}

static double update_nodes(nb_mshpoly_t *msh, uint32_t *N_adj,
			   uint32_t **adj, char *mask,
			   double (*density)(const double[2],
					     const void *data),
			   const void *density_data)
{
	uint32_t N_nodes = nb_mshpoly_get_N_nodes(msh);
	double max_d2 = 0;
	for (uint32_t i = 0; i < N_nodes; i++) {
		if (0 == mask[i]) {
			double p[2];
			get_new_node_position(msh, i, N_adj, adj, p,
					      density, density_data);
			double d2 = nb_utils2D_get_dist2(p, &(msh->nod[i*2]));
			if (d2 > max_d2)
				max_d2 = d2;
			msh->nod[i * 2] = p[0];
			msh->nod[i*2+1] = p[1];
		}
	}

	uint32_t N_sgm = nb_mshpoly_get_N_insgm(msh);
	for (uint32_t i = 0; i < N_sgm; i++) {
		uint32_t N_sgm_nodes = nb_mshpoly_insgm_get_N_nodes(msh, i);
		double sn[2];
		double s1[2];
		if (N_sgm_nodes > 0) {
			uint32_t id1 =
				nb_mshpoly_insgm_get_node(msh, i, 0);
			uint32_t id2 =
				nb_mshpoly_insgm_get_node(msh, i,
							  N_sgm_nodes - 1);
			s1[0] = nb_mshpoly_node_get_x(msh, id1);
			s1[1] = nb_mshpoly_node_get_y(msh, id1);
			double s2[2];
			s2[0] = nb_mshpoly_node_get_x(msh, id2);
			s2[1] = nb_mshpoly_node_get_y(msh, id2);
			double len = nb_utils2D_get_dist(s1, s2);
			sn[0] = (s2[0] - s1[0]) / len;
			sn[1] = (s2[1] - s1[1]) / len;			
		}
		for (uint32_t j = 0; j < N_sgm_nodes; j++) {
			uint32_t id = nb_mshpoly_insgm_get_node(msh, i, j);
			if (1 == mask[id]) {
				double p[2];
				get_new_node_position(msh, id, N_adj, adj, p,
						      density, density_data);
				double dot = (p[0] - s1[0]) * sn[0] +
					(p[1] - s1[1]) * sn[1];
				p[0] = s1[0] + dot * sn[0];
				p[1] = s1[1] + dot * sn[1];
				double *node = &(msh->nod[id*2]);
				double d2 = nb_utils2D_get_dist2(p, node);
				if (d2 > max_d2)
					max_d2 = d2;
				msh->nod[id * 2] = p[0];
				msh->nod[id*2+1] = p[1];
			}
		}
	}
	return max_d2;
}

static void get_new_node_position(nb_mshpoly_t *msh, uint32_t id,
				  uint32_t *N_adj, uint32_t **adj,
				  double new_p[2],
				  double (*density)(const double[2],
					     const void *data),
				  const void *density_data)
{
	double x = 0;
	double y = 0;
	for (int j = 0; j < N_adj[id]; j++) {
		uint32_t elem_id = adj[id][j];
		double den = 1;
		if (density != NULL)
			den = density(&(msh->cen[elem_id * 2]), density_data);
		x += msh->cen[elem_id * 2] * den;
		y += msh->cen[elem_id*2+1] * den;
	}
	new_p[0] = x / N_adj[id];
	new_p[1] = y / N_adj[id];
}
