#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>
#include <alloca.h>
#include <math.h>

#include "nb/math_bot.h"
#include "nb/graph_bot.h"
#include "nb/geometric_bot/utils2D.h"
#include "nb/geometric_bot/knn/bins2D.h"
#include "nb/geometric_bot/knn/bins2D_iterator.h"
#include "nb/geometric_bot/mesh/mesh2D.h"
#include "nb/geometric_bot/mesh/modules2D/graph_generator.h"
#include "nb/geometric_bot/mesh/partition/elements2D/mshquad.h"

#include "mshquad_struct.h"
#include "../../mesh2D_structs.h"

#define MAX(a,b) (((a)>(b))?(a):(b))

typedef struct {
	uint32_t N_matchs;
	uint32_t N_unmatched_trg;
} match_data;

static void set_arrays_memory(nb_mshquad_t *quad);
static void copy_nodes(nb_mshquad_t* quad, const nb_mshquad_t *const src_quad);
static void copy_edges(nb_mshquad_t* quad, const nb_mshquad_t *const src_quad);
static void copy_elems(nb_mshquad_t* quad, const nb_mshquad_t *const src_quad);
static void copy_vtx(nb_mshquad_t* quad, const nb_mshquad_t *const src_quad);
static void copy_N_nod_x_sgm(nb_mshquad_t* quad,
			     const nb_mshquad_t *const src_quad);
static uint32_t get_size_of_nod_x_sgm(const nb_mshquad_t *const quad);
static void set_mem_of_nod_x_sgm(nb_mshquad_t *quad, uint32_t memsize);
static void copy_nod_x_sgm(nb_mshquad_t* quad,
			   const nb_mshquad_t *const src_quad);
static void* malloc_quad(void);
static void set_quad_quality_as_weights(const nb_mesh_t *const mesh,
					nb_graph_t *graph);
static double get_quality(const msh_trg_t *trg1,
			  const msh_trg_t *trg2);
static bool adj_is_input_sgm(const msh_trg_t *trg1,
			     const msh_trg_t *trg2);
static void get_quad_matching_vtx(const msh_trg_t *const trg,
				  const msh_trg_t *const match_trg,
				  msh_vtx_t *quad_vtx[4]);
static double get_max_angle_distortion(msh_vtx_t* vtx[4]);
static void get_angle_vertices(msh_vtx_t* vtx[4],
			       double a[2], double b[2],
			       double c[2], uint32_t id);
static msh_trg_t * get_trg_adj(const msh_trg_t *const trg,
			       uint32_t id_adj);
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
static void init_elems(nb_mshquad_t *quad, const nb_mesh_t *const mesh,
		       const uint32_t *const matches,
		       uint32_t *new_elem_id);
static void set_trg_element(nb_mshquad_t *quad,
			    const msh_trg_t *const trg, uint32_t elem_id);
static void set_quad_element(nb_mshquad_t *quad,
			     const msh_trg_t *const trg, 
			     uint32_t match_id, uint32_t elem_id);
static msh_trg_t *get_match_trg(const msh_trg_t *const trg,
				uint32_t match_id);
static void set_quad_from_trg(nb_mshquad_t *quad,
			      const msh_trg_t *const trg,
			      const msh_trg_t *const match_trg,
			      uint32_t elem_id);
static void update_neighbors_ids(nb_mshquad_t *quad,
				 const uint32_t *const new_elem_id);
static void set_vtx(nb_mshquad_t *quad, const nb_mesh_t *const mesh);
static uint32_t set_N_nod_x_sgm(nb_mshquad_t *quad, const nb_mesh_t *const mesh);
static void set_nod_x_sgm(nb_mshquad_t *quad, const nb_mesh_t *const mesh);
static void set_sgm_nodes(nb_mshquad_t *quad,
			  const vcn_mesh_t *const mesh,
			  uint32_t sgm_id);
static void assemble_sgm_wire(nb_mshquad_t *quad, uint32_t sgm_id,
			      msh_edge_t *sgm_prev, msh_edge_t *sgm);

uint32_t nb_mshquad_get_memsize(void)
{
	return sizeof(nb_mshquad_t);
}

void nb_mshquad_init(void *mshquad_ptr)
{
	memset(mshquad_ptr, 0, nb_mshquad_get_memsize());
}

void nb_mshquad_copy(void *dest, const void *const src)
{
	memcpy(dest, src, nb_mshquad_get_memsize());
	nb_mshquad_t *quad = dest;
	const nb_mshquad_t *const src_quad = src;
	
	if (quad->N_elems > 0) {
		set_arrays_memory(quad);

		copy_nodes(quad, src_quad);
		copy_edges(quad, src_quad);
		copy_elems(quad, src_quad);
		copy_vtx(quad, src_quad);
		copy_N_nod_x_sgm(quad, src_quad);

		uint32_t memsize = get_size_of_nod_x_sgm(quad);
		set_mem_of_nod_x_sgm(quad, memsize);

		copy_nod_x_sgm(quad, src_quad);
	}
}

static void set_arrays_memory(nb_mshquad_t *quad)
{
	uint32_t nod_size = quad->N_nod * 2 * sizeof(*(quad->nod));
	uint32_t edg_size = quad->N_edg * 2 * sizeof(*(quad->edg));
	uint32_t type_size = quad->N_elems * sizeof(*(quad->type));
	uint32_t adj_size = quad->N_elems * 4 * sizeof(*(quad->adj));
	uint32_t ngb_size = quad->N_elems * 4 * sizeof(*(quad->ngb));
	uint32_t vtx_size = quad->N_vtx * sizeof(*(quad->vtx));
	uint32_t N_nod_x_sgm_size = quad->N_sgm *
		sizeof(*(quad->N_nod_x_sgm));
	uint32_t nod_x_sgm_size = quad->N_sgm * sizeof(*(quad->nod_x_sgm));

	uint32_t size = nod_size + edg_size + type_size + adj_size +
		ngb_size + vtx_size + N_nod_x_sgm_size + nod_x_sgm_size;

	char *memblock = malloc(size);

	quad->nod = (void*) memblock;
	quad->edg = (void*) ((char*)(quad->nod) + nod_size);
	quad->type = (void*) ((char*)(quad->edg) + edg_size);
	quad->adj = (void*) ((char*)(quad->type) + type_size);
	quad->ngb = (void*) ((char*)(quad->adj) + adj_size);
	quad->vtx = (void*) ((char*)(quad->ngb) + ngb_size);
	quad->N_nod_x_sgm = (void*) ((char*)(quad->vtx) + vtx_size);
	quad->nod_x_sgm = (void*) ((char*)(quad->N_nod_x_sgm) +
				   N_nod_x_sgm_size);
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

static void copy_N_nod_x_sgm(nb_mshquad_t* quad,
			     const nb_mshquad_t *const src_quad)
{
	memcpy(quad->N_nod_x_sgm, src_quad->N_nod_x_sgm,
	       quad->N_sgm * sizeof(*(quad->N_nod_x_sgm)));
}
static uint32_t get_size_of_nod_x_sgm(const nb_mshquad_t *const quad)
{
	uint32_t size = 0;
	for (uint32_t i = 0; i < quad->N_sgm; i++)
		size += quad->N_nod_x_sgm[i] *
			sizeof(**(quad->nod_x_sgm));
	return size;
}

static void set_mem_of_nod_x_sgm(nb_mshquad_t *quad, uint32_t memsize)
{
	char *memblock = malloc(memsize);
	for (uint32_t i = 0; i < quad->N_sgm; i++) {
		quad->nod_x_sgm[i] = (void*) memblock;
		memblock += quad->N_nod_x_sgm[i] *
			sizeof(**(quad->nod_x_sgm));
	}
}

static void copy_nod_x_sgm(nb_mshquad_t* quad,
			   const nb_mshquad_t *const src_quad)
{	
	for (int i = 0; i < quad->N_sgm; i++) {
		memcpy(&(quad->nod_x_sgm[i]), &(src_quad->nod_x_sgm[i]),
		       quad->N_nod_x_sgm[i] * sizeof(**(quad->nod_x_sgm)));
	}

}

void nb_mshquad_finish(void *mshquad_ptr)
{
	nb_mshquad_clear(mshquad_ptr);
}

void* nb_mshquad_create(void)
{
	nb_mshquad_t *quad = malloc_quad();
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
	nb_mshquad_t *quad = malloc_quad();
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
		free(quad->nod_x_sgm[0]);
		free(quad->nod);		
	}
	memset(mshquad_ptr, 0, nb_mshquad_get_memsize());	
}

uint32_t nb_mshquad_get_N_invtx(const void *msh)
{
	const nb_mshquad_t *mshquad = msh;
	return mshquad->N_vtx;
}

uint32_t nb_mshquad_get_N_insgm(const void *msh)
{
	const nb_mshquad_t *mshquad = msh;
	return mshquad->N_sgm;
}

uint32_t nb_mshquad_get_N_nodes(const void *msh)
{
	const nb_mshquad_t *mshquad = msh;
	return mshquad->N_nod;
}

uint32_t nb_mshquad_get_N_edges(const void *msh)
{
	const nb_mshquad_t *mshquad = msh;
	return mshquad->N_edg;
}

uint32_t nb_mshquad_get_N_elems(const void *msh)
{
	const nb_mshquad_t *mshquad = msh;
	return mshquad->N_elems;
}

double nb_mshquad_node_get_x(const void *msh, uint32_t id)
{
	const nb_mshquad_t *mshquad = msh;
	return mshquad->nod[id * 2];
}

double nb_mshquad_node_get_y(const void *msh, uint32_t id)
{
	const nb_mshquad_t *mshquad = msh;
	return mshquad->nod[id*2+1];
}

uint32_t nb_mshquad_edge_get_1n(const void *msh, uint32_t id)
{
	const nb_mshquad_t *mshquad = msh;
	return mshquad->edg[id * 2];
}

uint32_t nb_mshquad_edge_get_2n(const void *msh, uint32_t id)
{
	const nb_mshquad_t *mshquad = msh;
	return mshquad->edg[id*2+1];
}

void nb_mshquad_edge_get_midpoint(const void *msh,
				  uint32_t face_id,
				  double w, double midpoint[2])
{
	const nb_mshquad_t *quad = msh;

	uint32_t n1 =  quad->edg[face_id * 2];
	uint32_t n2 =  quad->edg[face_id*2+1];

	double *s1 = &(quad->nod[n1 * 2]);
	double *s2 = &(quad->nod[n2 * 2]);

	midpoint[0] = (1 - w) * s1[0] + w * s2[0];
	midpoint[1] = (1 - w) * s1[1] + w * s2[1];
}

double nb_mshquad_edge_get_normal(const void *msh, uint32_t face_id,
				  double normal[2])
{

	const nb_mshquad_t *quad = msh;
	uint32_t n1 =  quad->edg[face_id * 2];
	uint32_t n2 =  quad->edg[face_id*2+1];
	double *s1 = &(quad->nod[n1 * 2]);
	double *s2 = &(quad->nod[n2 * 2]);
	double length = vcn_utils2D_get_dist(s1, s2);
	normal[0] =  (s2[1] - s1[1]) / length;
	normal[1] = -(s2[0] - s1[0]) / length;
	return length;
}

double nb_mshquad_elem_get_x(const void *msh, uint32_t id)
{
	uint32_t idx = nb_mshquad_elem_get_adj(msh, id, 0);
	double x = nb_mshquad_node_get_x(msh, idx);

	idx = nb_mshquad_elem_get_adj(msh, id, 1);
	x += nb_mshquad_node_get_x(msh, idx);

	idx = nb_mshquad_elem_get_adj(msh, id, 2);
	x += nb_mshquad_node_get_x(msh, idx);

	double div = 3.0;
	if (nb_mshquad_elem_is_quad(msh, id)) {
		div = 4.0;
		idx = nb_mshquad_elem_get_adj(msh, id, 3);
		x += nb_mshquad_node_get_x(msh, idx);
	}
	return x / div;
}

double nb_mshquad_elem_get_y(const void *msh, uint32_t id)
{
	uint32_t idx = nb_mshquad_elem_get_adj(msh, id, 0);
	double y = nb_mshquad_node_get_y(msh, idx);

	idx = nb_mshquad_elem_get_adj(msh, id, 1);
	y += nb_mshquad_node_get_y(msh, idx);

	idx = nb_mshquad_elem_get_adj(msh, id, 2);
	y += nb_mshquad_node_get_y(msh, idx);

	double div = 3.0;
	if (nb_mshquad_elem_is_quad(msh, id)) {
		div = 4.0;
		idx = nb_mshquad_elem_get_adj(msh, id, 3);
		y += nb_mshquad_node_get_y(msh, idx);
	}
	return y / div;
}

double nb_mshquad_elem_get_area(const void *msh, uint32_t id)
{
	const nb_mshquad_t *mshquad = msh;

	uint32_t v1 = nb_mshquad_elem_get_adj(msh, id, 0);
	uint32_t v2 = nb_mshquad_elem_get_adj(msh, id, 1);
	uint32_t v3 = nb_mshquad_elem_get_adj(msh, id, 2);

	double *t1 = &(mshquad->nod[v1 * 2]);
	double *t2 = &(mshquad->nod[v2 * 2]);
	double *t3 = &(mshquad->nod[v3 * 2]);

	double area = vcn_utils2D_get_trg_area(t1, t2, t3);

	if (nb_mshquad_elem_is_quad(msh, id)) {
		uint32_t v4 = nb_mshquad_elem_get_adj(msh, id, 3);
		double *t4 = &(mshquad->nod[v4 * 2]);
		area += vcn_utils2D_get_trg_area(t1, t3, t4);
	}

	return area;
}

double nb_mshquad_elem_get_radius(const void *msh, uint32_t id)
{
	const nb_mshquad_t *quad = msh;
	double x[2];
	x[0] = nb_mshquad_elem_get_x(msh, id);
	x[1] = nb_mshquad_elem_get_y(msh, id);
	uint32_t N_adj = nb_mshquad_elem_get_N_adj(msh, id);
	double max = 0;
	for (uint16_t i = 0; i < N_adj; i++) {
		uint32_t nid = nb_mshquad_elem_get_adj(msh, id, i);
		double *ni = &(quad->nod[nid * 2]);
		double dist2 = vcn_utils2D_get_dist2(x, ni);
		if (dist2 > max)
			max = dist2;
	}
	return sqrt(max);	
}

double nb_mshquad_elem_get_apotem(const void *msh, uint32_t id)
{
	const nb_mshquad_t *quad = msh;
	double x[2];
	x[0] = nb_mshquad_elem_get_x(msh, id);
	x[1] = nb_mshquad_elem_get_y(msh, id);
	uint32_t N_adj = nb_mshquad_elem_get_N_adj(msh, id);
	double max = 0;
	for (uint16_t i = 0; i < N_adj; i++) {
		uint32_t id1 = nb_mshquad_elem_get_adj(msh, id, i);
		uint32_t id2 = nb_mshquad_elem_get_adj(msh, id, (i+1)%N_adj);
		double *n1 = &(quad->nod[id1 * 2]);
		double *n2 = &(quad->nod[id2 * 2]);
		double closest[2];
		vcn_utils2D_get_closest_pnt_to_sgm(n1, n2, x, closest);
		double dist2 = vcn_utils2D_get_dist2(x, closest);
		if (dist2 > max)
			max = dist2;
	}
	return sqrt(max);
}

uint32_t nb_mshquad_elem_find_edge(const void *msh, uint32_t id,
				   uint16_t local_face_id)
{
	const nb_mshquad_t *quad = msh;
	uint16_t N_adj = nb_mshquad_elem_get_N_adj(msh, id);

	uint16_t l1 = local_face_id;
	uint16_t l2 = (local_face_id + 1) % N_adj;

	uint32_t n1 = quad->adj[id*4+l1];
	uint32_t n2 = quad->adj[id*4+l2];

	uint32_t N_edges = quad->N_edg;
	uint32_t edge_id = N_edges;
	for (uint32_t i = 0; i < N_edges; i++) {
		uint32_t v1 = quad->edg[i * 2];
		uint32_t v2 = quad->edg[i*2+1];
		if ((n1 == v1 && n2 == v2) || (n1 == v2 && n2 == v1)) {
			edge_id = i;
			break;
		}
	}
	return edge_id;
}

double nb_mshquad_elem_face_get_length(const void *msh,
				       uint32_t elem_id,
				       uint16_t face_id)
{
	const nb_mshquad_t *mshquad = msh;

	uint16_t N_adj = nb_mshquad_elem_get_N_adj(msh, elem_id);
	uint32_t v1 = nb_mshquad_elem_get_adj(msh, elem_id, face_id);
	uint32_t v2 = nb_mshquad_elem_get_adj(msh, elem_id,
					      (face_id + 1) % N_adj);

	double *n1 = &(mshquad->nod[v1 * 2]);
	double *n2 = &(mshquad->nod[v2 * 2]);

	return vcn_utils2D_get_dist(n1, n2);	
}

double nb_mshquad_elem_face_get_normal(const void *msh, uint32_t elem_id,
				       uint16_t face_id, double normal[2])
{
	const nb_mshquad_t *mshquad = msh;
	uint16_t N_adj = nb_mshquad_elem_get_N_adj(msh, elem_id);
	uint32_t n1 = nb_mshquad_elem_get_adj(msh, elem_id, face_id);
	uint32_t n2 = nb_mshquad_elem_get_adj(msh, elem_id,
					      (face_id + 1) % N_adj);
	double *s1 = &(mshquad->nod[n1 * 2]);
	double *s2 = &(mshquad->nod[n2 * 2]);
	double length = vcn_utils2D_get_dist(s1, s2);
	normal[0] =  (s2[1] - s1[1]) / length;
	normal[1] = -(s2[0] - s1[0]) / length;
	return length;
}

double nb_mshquad_elem_ngb_get_normal(const void *msh, uint32_t elem_id,
				      uint16_t ngb_id, double normal[2])
{
	const nb_mshquad_t *mshquad = msh;
	uint32_t N_elems = nb_mshquad_get_N_elems(msh);
	uint32_t nid = nb_mshquad_elem_get_ngb(msh, elem_id, ngb_id);
	memset(normal, 0, 2 * sizeof(double));
	double dist = 0;
	if (nid < N_elems) {
		double id1[2];
		id1[0] = nb_mshquad_elem_get_x(msh, elem_id);
		id1[1] = nb_mshquad_elem_get_y(msh, elem_id);

		double id2[2];
		id2[0] = nb_mshquad_elem_get_x(msh, nid);
		id2[1] = nb_mshquad_elem_get_y(msh, nid);

		dist = vcn_utils2D_get_dist(id1, id2);
		normal[0] = (id2[0] - id1[0]) / dist;
		normal[1] = (id2[1] - id1[1]) / dist;
	}
	return dist;
}

uint32_t nb_mshquad_elem_get_N_adj(const void *msh, uint32_t id)
{
	uint32_t out;
	if (nb_mshquad_elem_is_quad(msh, id))
		out = 4;
	else
		out = 3;
	return out;
}

uint32_t nb_mshquad_elem_get_adj(const void *msh, uint32_t elem_id,
				 uint8_t adj_id)
{
	const nb_mshquad_t *mshquad = msh;
	return mshquad->adj[elem_id * 4 + adj_id];
}

uint32_t nb_mshquad_elem_get_ngb(const void *msh,
				 uint32_t elem_id, uint8_t ngb_id)
{
	const nb_mshquad_t *mshquad = msh;
	return mshquad->ngb[elem_id * 4 + ngb_id];
}

bool nb_mshquad_elem_has_ngb(const void *msh, uint32_t elem_id,
			     uint16_t ngb_id)
{
	uint32_t N_elems = nb_mshquad_get_N_elems(msh);
	uint32_t id = nb_mshquad_elem_get_ngb(msh, elem_id, ngb_id);
	return id < N_elems;
}

bool nb_mshquad_elem_is_quad(const void *msh, uint32_t elem_id)
{
	const nb_mshquad_t *mshquad = msh;
	return 0 == mshquad->type[elem_id];
}

uint32_t nb_mshquad_get_invtx(const void *msh, uint32_t id)
{
	const nb_mshquad_t *mshquad = msh;
	return mshquad->vtx[id];
}

uint32_t nb_mshquad_insgm_get_N_nodes(const void *msh, uint32_t id)
{
	const nb_mshquad_t *mshquad = msh;
	return mshquad->N_nod_x_sgm[id];
}

uint32_t nb_mshquad_insgm_get_node(const void *msh, uint32_t sgm_id,
				     uint32_t node_id)
{
	const nb_mshquad_t *mshquad = msh;
	return mshquad->nod_x_sgm[sgm_id][node_id];
}

void nb_mshquad_load_from_mesh(void *mshquad, nb_mesh_t *mesh)
{
	if (vcn_mesh_get_N_trg(mesh) > 0) {
		mesh_enumerate_vtx((vcn_mesh_t*)mesh);
		mesh_enumerate_trg((vcn_mesh_t*)mesh);
		nb_graph_t *graph = vcn_mesh_create_elem_graph(mesh);
		nb_graph_init_edge_weights(graph);

		set_quad_quality_as_weights(mesh, graph);

		uint32_t *matches = malloc(graph->N * sizeof(*matches));
		nb_graph_matching_greedy(graph, matches);
	
		set_mshquad(mshquad, graph, mesh, matches);

		free(matches);

		nb_graph_finish_edge_weights(graph);
		nb_graph_destroy(graph);
	}
}

static void set_quad_quality_as_weights(const nb_mesh_t *const mesh,
					nb_graph_t *graph)
{
	uint16_t iter_size = nb_iterator_get_memsize();
	nb_iterator_t *iter = alloca(iter_size);
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, mesh->ht_trg);
	uint32_t elem_id = 0;
	while (nb_iterator_has_more(iter)) {
		msh_trg_t* trg = (msh_trg_t*) nb_iterator_get_next(iter);
		uint32_t id = trg->id;

		for (uint32_t i = 0; i < graph->N_adj[id]; i++) {
			uint32_t id_adj = graph->adj[id][i];
			msh_trg_t *trg_adj = get_trg_adj(trg, id_adj);
			graph->wij[id][i] = get_quality(trg, trg_adj);
		}
	}
	nb_iterator_finish(iter);	
}

static double get_quality(const msh_trg_t *trg1,
			  const msh_trg_t *trg2)
{
	double quality;
	if (adj_is_input_sgm(trg1, trg2)) {
		quality = 0.0;
	} else {
		msh_vtx_t *quad_vtx[4];
		get_quad_matching_vtx(trg1, trg2, quad_vtx);
		double maxk = get_max_angle_distortion(quad_vtx);
		quality = MAX(0.0, 1.0 - 2.0 * maxk / NB_PI);
	}
	return quality;
}

static bool adj_is_input_sgm(const msh_trg_t *trg1,
			     const msh_trg_t *trg2)
{
	bool out;
	if (trg1->t1 == trg2) {
		msh_edge_t *edg = trg1->s1;
		out = medge_is_subsgm(edg);
	} else if (trg1->t2 == trg2) {
		msh_edge_t *edg = trg1->s2;
		out = medge_is_subsgm(edg);
	} else if (trg1->t3 == trg2) {
		msh_edge_t *edg = trg1->s3;
		out = medge_is_subsgm(edg);
	} else {
		out = false;
	}
	return out;
}

static void get_quad_matching_vtx(const msh_trg_t *const trg,
				  const msh_trg_t *const match_trg,
				  msh_vtx_t *quad_vtx[4])
{
	if (trg->t1 == match_trg) {
		quad_vtx[0] = trg->v1;
		quad_vtx[1] = mtrg_get_opposite_vertex(match_trg, trg->s1);
		quad_vtx[2] = trg->v2;
		quad_vtx[3] = trg->v3;
	} else if (trg->t2 == match_trg) {
		quad_vtx[0] = trg->v2;
		quad_vtx[1] = mtrg_get_opposite_vertex(match_trg, trg->s2);
		quad_vtx[2] = trg->v3;
		quad_vtx[3] = trg->v1;
	} else {
		quad_vtx[0] = trg->v3;
		quad_vtx[1] = mtrg_get_opposite_vertex(match_trg, trg->s3);
		quad_vtx[2] = trg->v1;
		quad_vtx[3] = trg->v2;
	}
}

static double get_max_angle_distortion(msh_vtx_t* vtx[4])
{
	double max_distortion = 0.0;
	for (uint32_t i = 0; i < 4; i++) {
		double a[2], b[2], c[2];
		get_angle_vertices(vtx, a, b, c, i);
		double angle = nb_utils2D_get_2vec_angle(a, b, c);
		double distortion = fabs(NB_PI/2.0 - angle);
		max_distortion = MAX(max_distortion, distortion);
	}
	return max_distortion;
}

static void get_angle_vertices(msh_vtx_t* vtx[4],
			       double a[2], double b[2],
			       double c[2], uint32_t id)
{
	uint16_t vtx_size = 2 * sizeof(double);
	if (0 == id) {
		memcpy(a, vtx[3]->x, vtx_size);
		memcpy(b, vtx[0]->x, vtx_size);
		memcpy(c, vtx[1]->x, vtx_size);
	} else if (1 == id) {
		memcpy(a, vtx[0]->x, vtx_size);
		memcpy(b, vtx[1]->x, vtx_size);
		memcpy(c, vtx[2]->x, vtx_size);
	} else if (2 == id) {
		memcpy(a, vtx[1]->x, vtx_size);
		memcpy(b, vtx[2]->x, vtx_size);
		memcpy(c, vtx[3]->x, vtx_size);
	} else {
		memcpy(a, vtx[2]->x, vtx_size);
		memcpy(b, vtx[3]->x, vtx_size);
		memcpy(c, vtx[0]->x, vtx_size);
	}
}

static msh_trg_t* get_trg_adj(const msh_trg_t *const trg,
			      uint32_t id_adj)
{
	msh_trg_t *trg_adj = NULL;
	if (NULL != trg->t1) {
		if (id_adj == trg->t1->id)
			trg_adj = trg->t1;
	}
	if (NULL == trg_adj && NULL != trg->t2) {
		if (id_adj == trg->t2->id)
			trg_adj = trg->t2;
	}
	if (NULL == trg_adj && NULL != trg->t3) {
		if (id_adj == trg->t3->id)
			trg_adj = trg->t3;
	}
	return trg_adj;
}

static void set_mshquad(nb_mshquad_t *quad,
			const nb_graph_t *const graph,
			const nb_mesh_t *const mesh,
			const uint32_t *const matches)
{
	match_data *data = alloca(sizeof(match_data));
	get_match_data(graph, matches, data);
	
	quad->N_nod = vcn_mesh_get_N_vtx(mesh);
	quad->N_edg = vcn_mesh_get_N_edg(mesh) - data->N_matchs;
	quad->N_elems = data->N_matchs + data->N_unmatched_trg;
	quad->N_vtx = mesh->N_input_vtx;
	quad->N_sgm = mesh->N_input_sgm;
	
	set_arrays_memory(quad);

	set_nodes(quad, mesh);
	set_edges(quad, mesh, matches);
	set_elems(quad, mesh, matches);
	set_vtx(quad, mesh);

	uint32_t nod_x_sgm_memsize = set_N_nod_x_sgm(quad, mesh);
	set_mem_of_nod_x_sgm(quad, nod_x_sgm_memsize);
	set_nod_x_sgm(quad, mesh);
}

static void set_nodes(nb_mshquad_t *quad, const nb_mesh_t *const mesh)
{
	vcn_bins2D_iter_t* iter = alloca(vcn_bins2D_iter_get_memsize());
	vcn_bins2D_iter_init(iter);
	vcn_bins2D_iter_set_bins(iter, mesh->ug_vtx);
	while (vcn_bins2D_iter_has_more(iter)) {
		const msh_vtx_t* vtx = vcn_bins2D_iter_get_next(iter);
		uint32_t id = mvtx_get_id(vtx);
		quad->nod[id * 2] = vtx->x[0] / mesh->scale + mesh->xdisp;
		quad->nod[id*2+1] = vtx->x[1] / mesh->scale + mesh->ydisp;
	}
	vcn_bins2D_iter_finish(iter);	
}

static void set_edges(nb_mshquad_t *quad, const nb_mesh_t *const mesh,
		      const uint32_t *const matches)
{
	uint32_t i = 0;
	nb_iterator_t *iter = alloca(nb_iterator_get_memsize());
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, mesh->ht_edge);
	while (nb_iterator_has_more(iter)) {
		msh_edge_t *edge = (msh_edge_t*) nb_iterator_get_next(iter);
		if (edge_is_not_matched(edge, matches)) {
			quad->edg[i * 2] = mvtx_get_id(edge->v1);
			quad->edg[i*2+1] = mvtx_get_id(edge->v2);
			i += 1;
		}
	}
	nb_iterator_finish(iter);
}

static bool edge_is_not_matched(const msh_edge_t *const edge,
				const uint32_t *const matches)
{
	bool out = true;
	if (!medge_is_boundary(edge)) {
		uint32_t id1 = edge->t1->id;
		uint32_t id2 = edge->t2->id;
		if (matches[id1] == id2)
			out = false;
	}
	return out;
}

static void set_elems(nb_mshquad_t *quad, const nb_mesh_t *const mesh,
		      const uint32_t *const matches)
{
	uint32_t N_trg = vcn_mesh_get_N_trg(mesh);
	uint32_t memsize = N_trg * sizeof(uint32_t);
	uint32_t *new_elem_id = NB_SOFT_MALLOC(memsize);
	init_elems(quad, mesh, matches, new_elem_id);
	update_neighbors_ids(quad, new_elem_id);
	NB_SOFT_FREE(memsize, new_elem_id);
}

static void init_elems(nb_mshquad_t *quad, const nb_mesh_t *const mesh,
		       const uint32_t *const matches,
		       uint32_t *new_elem_id)
{
	uint16_t iter_size = nb_iterator_get_memsize();
	nb_iterator_t *iter = alloca(iter_size);
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, mesh->ht_trg);
	uint32_t elem_id = 0;
	while (nb_iterator_has_more(iter)) {
		msh_trg_t* trg = (msh_trg_t*) nb_iterator_get_next(iter);
		uint32_t id = trg->id;
		uint32_t match_id = matches[id];

		if (match_id == id) {
			set_trg_element(quad, trg, elem_id);
			new_elem_id[id] = elem_id;
			elem_id += 1;
		} else {
			if (id < match_id) {
				set_quad_element(quad, trg,
						 match_id, elem_id);
				new_elem_id[id] = elem_id;
				new_elem_id[match_id] = elem_id;
				elem_id += 1;
			}
		}
	}
	nb_iterator_finish(iter);
}

static void set_trg_element(nb_mshquad_t *quad,
			    const msh_trg_t *const trg, uint32_t elem_id)
{
	quad->type[elem_id] = 1;

	quad->adj[elem_id * 4] = mvtx_get_id(trg->v1);
	quad->adj[elem_id*4+1] = mvtx_get_id(trg->v2);
	quad->adj[elem_id*4+2] = mvtx_get_id(trg->v3);
	quad->adj[elem_id*4+3] = quad->N_nod;

	uint32_t t1;
	if (NULL != trg->t1)
		t1 = trg->t1->id;
	else
		t1 = UINT32_MAX;

	uint32_t t2;
	if (NULL != trg->t2)
		t2 = trg->t2->id;
	else
		t2 = UINT32_MAX;

	uint32_t t3;
	if (NULL != trg->t3)
		t3 = trg->t3->id;
	else
		t3 = UINT32_MAX;

	quad->ngb[elem_id * 4] = t1;
	quad->ngb[elem_id*4+1] = t2;
	quad->ngb[elem_id*4+2] = t3;
	quad->ngb[elem_id*4+3] = UINT32_MAX;
}

static void set_quad_element(nb_mshquad_t *quad,
			     const msh_trg_t *const trg, 
			     uint32_t match_id, uint32_t elem_id)
{
	quad->type[elem_id] = 0;

	msh_trg_t *match_trg = get_match_trg(trg, match_id);
	set_quad_from_trg(quad, trg, match_trg, elem_id);
}

static msh_trg_t *get_match_trg(const msh_trg_t *const trg,
				uint32_t match_id)
{
	msh_trg_t *match_trg = NULL;
	if (NULL != trg->t1) {
		if (match_id == trg->t1->id)
			match_trg = trg->t1;
	}
	if (NULL == match_trg && NULL != trg->t2) {
		if (match_id == trg->t2->id)
			match_trg = trg->t2;
	}
	if (NULL == match_trg && NULL != trg->t3) {
		if (match_id == trg->t3->id)
			match_trg = trg->t3;
	}
	return match_trg;
}

static void set_quad_from_trg(nb_mshquad_t *quad,
			      const msh_trg_t *const trg,
			      const msh_trg_t *const match_trg,
			      uint32_t elem_id)
{
	msh_vtx_t *vtx[4];
	get_quad_matching_vtx(trg, match_trg, vtx);
	
	msh_trg_t *trg1 = mtrg_get_left_triangle(match_trg, vtx[1]);
	msh_trg_t *trg2 = mtrg_get_right_triangle(match_trg, vtx[1]);
	msh_trg_t *trg3, *trg4;
	if (trg->t1 == match_trg) {
		trg3 = trg->t2;
		trg4 = trg->t3;
	} else if (trg->t2 == match_trg) {
		trg3 = trg->t3;
		trg4 = trg->t1;
	} else {
		trg3 = trg->t1;
		trg4 = trg->t2;
	}
	quad->adj[elem_id * 4] = mvtx_get_id(vtx[0]);
	quad->adj[elem_id*4+1] = mvtx_get_id(vtx[1]);
	quad->adj[elem_id*4+2] = mvtx_get_id(vtx[2]);
	quad->adj[elem_id*4+3] = mvtx_get_id(vtx[3]);

	uint32_t t1;
	if (NULL != trg1)
		t1 = trg1->id;
	else
	        t1 = UINT32_MAX;

	uint32_t t2;
	if (NULL != trg2)
		t2 = trg2->id;
	else
	        t2 = UINT32_MAX;

	uint32_t t3;
	if (NULL != trg3)
		t3 = trg3->id;
	else
	        t3 = UINT32_MAX;

	uint32_t t4;
	if (NULL != trg4)
		t4 = trg4->id;
	else
	        t4 = UINT32_MAX;

	quad->ngb[elem_id * 4] = t1;
	quad->ngb[elem_id*4+1] = t2;
	quad->ngb[elem_id*4+2] = t3;
	quad->ngb[elem_id*4+3] = t4;
}

static void update_neighbors_ids(nb_mshquad_t *quad,
				 const uint32_t *const new_elem_id)
{
	for (uint32_t i = 0; i < 4 * quad->N_elems; i++) {
		if (quad->ngb[i] == UINT32_MAX)
			quad->ngb[i] = quad->N_elems;
		else
			quad->ngb[i] = new_elem_id[quad->ngb[i]];
	}
}

static void set_vtx(nb_mshquad_t *quad, const nb_mesh_t *const mesh)
{
	for (uint32_t i = 0; i < quad->N_vtx; i++) {
		if (NULL != mesh->input_vtx[i]) {
			msh_vtx_t *vtx = mesh->input_vtx[i];
			quad->vtx[i] = mvtx_get_id(vtx);
		} else {
			quad->vtx[i] = quad->N_nod;
		}
	}
}

static uint32_t set_N_nod_x_sgm(nb_mshquad_t *quad, const nb_mesh_t *const mesh)
{
	uint32_t N_nod = 0;
	for (uint32_t i = 0; i < quad->N_sgm; i++) {
		if (NULL != mesh->input_sgm[i]) {
			msh_edge_t* sgm = mesh->input_sgm[i];
			uint32_t counter = 0;
			while (NULL != sgm) {
				counter++;
				sgm = medge_subsgm_next(sgm);
			}
			quad->N_nod_x_sgm[i] = counter + 1;
			N_nod += quad->N_nod_x_sgm[i];
		} else {
			quad->N_nod_x_sgm[i] = 0;
		}
	}
	return N_nod * sizeof(**(quad->nod_x_sgm));
}

static void set_nod_x_sgm(nb_mshquad_t *quad, const nb_mesh_t *const mesh)
{
	for (uint32_t i = 0; i < quad->N_sgm; i++) {
		if (NULL != mesh->input_sgm[i])
			set_sgm_nodes(quad, mesh, i);
	}
}

static void set_sgm_nodes(nb_mshquad_t *quad,
			  const vcn_mesh_t *const mesh,
			  uint32_t sgm_id)
{
	msh_edge_t *sgm_prev = mesh->input_sgm[sgm_id];
	msh_edge_t *sgm = medge_subsgm_next(sgm_prev);
	if (NULL == sgm) {
		quad->nod_x_sgm[sgm_id][0] = mvtx_get_id(sgm_prev->v1);
		quad->nod_x_sgm[sgm_id][1] = mvtx_get_id(sgm_prev->v2);
	} else {
		assemble_sgm_wire(quad, sgm_id, sgm_prev, sgm);
	}
}

static void assemble_sgm_wire(nb_mshquad_t *quad, uint32_t sgm_id,
			      msh_edge_t *sgm_prev, msh_edge_t *sgm)
{
	uint32_t idx = 0;
	uint32_t id_chain;
	uint32_t id1 = mvtx_get_id(sgm_prev->v1);
	uint32_t id2 = mvtx_get_id(sgm_prev->v2);
	uint32_t id1n = mvtx_get_id(sgm->v1);
	uint32_t id2n = mvtx_get_id(sgm->v2);
	if (id2 == id1n || id2 == id2n) {
		quad->nod_x_sgm[sgm_id][idx++] =  id1;
		quad->nod_x_sgm[sgm_id][idx++] =  id2;
		id_chain = id2;
	} else {
		quad->nod_x_sgm[sgm_id][idx++] =  id2;
		quad->nod_x_sgm[sgm_id][idx++] =  id1;
		id_chain = id1;
	}
	while (NULL != sgm) {
		sgm_prev = sgm;
		uint32_t id1 = mvtx_get_id(sgm_prev->v1);
		uint32_t id2 = mvtx_get_id(sgm_prev->v2);
		if (id1 == id_chain) {
			quad->nod_x_sgm[sgm_id][idx++] =  id2;
			id_chain = id2;
		} else {
			quad->nod_x_sgm[sgm_id][idx++] =  id1;
			id_chain = id1;
		}
		sgm = medge_subsgm_next(sgm);
	}
}

static void get_match_data(const nb_graph_t *const graph,
			   const uint32_t *const matches,
			   match_data *data)
{
	data->N_matchs = 0;
	data->N_unmatched_trg = 0;
	for (uint32_t i = 0; i < graph->N; i++) {
		if (matches[i] != i)
			data->N_matchs += 1;
		else
			data->N_unmatched_trg += 1;
	}
	data->N_matchs /= 2; /* Double counting */
}

void nb_mshquad_get_enveloping_box(const void *mshquad_ptr, double box[4])
{
	const nb_mshquad_t *mshquad = mshquad_ptr;
	vcn_utils2D_get_enveloping_box(mshquad->N_nod, mshquad->nod,
				       2 * sizeof(*(mshquad->nod)),
				       vcn_utils2D_get_x_from_darray,
				       vcn_utils2D_get_y_from_darray,
				       box);
}

bool nb_mshquad_is_vtx_inside(const void *msh, double x, double y)
{
	return false;/* PENDING */
}

void nb_mshquad_build_model(const void *msh, nb_model_t *model)
{
	;/* PENDING */
}

void nb_mshquad_build_model_disabled_elems(const void *msh,
					   const bool *elems_enabled,
					   nb_model_t *model,
					   uint32_t *N_input_vtx,
					   uint32_t **input_vtx)
{
	;/* PENDING */
}
