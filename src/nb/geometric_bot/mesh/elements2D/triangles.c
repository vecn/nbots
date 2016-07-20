#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <alloca.h>

#include "nb/container_bot/container.h"
#include "nb/container_bot/iterator.h"
#include "nb/geometric_bot/utils2D.h"
#include "nb/geometric_bot/knn/bins2D.h"
#include "nb/geometric_bot/knn/bins2D_iterator.h"
#include "nb/geometric_bot/mesh/elements2D/triangles_struct.h"
#include "nb/geometric_bot/mesh/elements2D/triangles.h"
#include "nb/geometric_bot/mesh/elements2D/trg_exporter.h"

#include "../mesh2D_structs.h"

static void msh3trg_set_N_vtx(void *msh3trg_ptr, uint32_t N);
static void msh3trg_malloc_vtx(void *msh3trg_ptr);
static void msh3trg_set_vtx(void *msh3trg_ptr, uint32_t i,
			    double x, double y);
static void msh3trg_set_N_edg(void *msh3trg_ptr, uint32_t N);
static void msh3trg_malloc_edg(void *msh3trg_ptr);
static void msh3trg_set_edg(void *msh3trg_ptr, uint32_t i,
			    uint32_t v1, uint32_t v2);
static void msh3trg_set_N_trg(void *msh3trg_ptr, uint32_t N);
static void msh3trg_malloc_trg(void *msh3trg_ptr, bool include_neighbours);
static void msh3trg_set_trg(void *msh3trg_ptr, uint32_t i,
			    uint32_t v1, uint32_t v2, uint32_t v3);
static void msh3trg_set_trg_neighbours(void *msh3trg_ptr,
				       uint32_t i, uint32_t t1,
				       uint32_t t2, uint32_t t3);
static void msh3trg_set_N_input_vtx(void *msh3trg_ptr, uint32_t N);
static void msh3trg_malloc_input_vtx(void *msh3trg_ptr);
static void msh3trg_set_input_vtx(void *msh3trg_ptr, uint32_t i,
				  uint32_t vtx_id);
static void msh3trg_set_N_input_sgm(void *msh3trg_ptr, uint32_t N);
static void msh3trg_malloc_input_sgm_table(void *msh3trg_ptr);
static void msh3trg_input_sgm_set_N_vtx(void *msh3trg_ptr, uint32_t i,
					uint32_t N);
static void msh3trg_input_sgm_malloc_vtx(void *msh3trg_ptr, uint32_t i);
static void msh3trg_input_sgm_set_vtx(void *msh3trg_ptr, uint32_t isgm,
				      uint32_t ivtx, uint32_t vtx_id);
static uint32_t itrg_get_right_triangle
                         (const vcn_msh3trg_t *const delaunay, 
			  const bool *const enabled_elements,
			  uint32_t itrg, uint32_t ivtx);
static uint32_t itrg_get_left_triangle
                         (const vcn_msh3trg_t *const delaunay,
			  const bool *const enabled_elements,
			  uint32_t itrg, uint32_t ivtx);


uint32_t vcn_msh3trg_get_memsize(void)
{
	return sizeof(vcn_msh3trg_t);
}

void vcn_msh3trg_init(vcn_msh3trg_t *msh3trg)
{
	uint32_t memsize = vcn_msh3trg_get_memsize();
	memset(msh3trg, 0, memsize);
}

void vcn_msh3trg_finish(vcn_msh3trg_t *msh3trg)
{
	vcn_msh3trg_clear(msh3trg);
}

vcn_msh3trg_t* vcn_msh3trg_create(void)
{
	uint32_t memsize = vcn_msh3trg_get_memsize();
	vcn_msh3trg_t* msh3trg = malloc(memsize);
	vcn_msh3trg_init(msh3trg);
	return msh3trg;
}

vcn_msh3trg_t* vcn_msh3trg_clone(vcn_msh3trg_t* msh3trg)
{
	vcn_msh3trg_t* clone = vcn_msh3trg_create();
	clone->N_vertices = msh3trg->N_vertices;
	clone->vertices = malloc(2 * clone->N_vertices * sizeof(*(clone->vertices)));
	memcpy(clone->vertices, msh3trg->vertices, 
	       2 * clone->N_vertices * sizeof(*(clone->vertices)));

	clone->N_triangles = msh3trg->N_triangles;
	if (clone->N_triangles > 0) {
		clone->vertices_forming_triangles =
			malloc(3 * clone->N_triangles *
			       sizeof(*(clone->vertices_forming_triangles)));
		memcpy(clone->vertices_forming_triangles,
		       msh3trg->vertices_forming_triangles,
		       3 * clone->N_triangles * 
		       sizeof(*(clone->vertices_forming_triangles)));
	}
	if (NULL != msh3trg->triangles_sharing_sides) {
		clone->triangles_sharing_sides =
			malloc(3 * clone->N_triangles *
			       sizeof(*(clone->triangles_sharing_sides)));
		memcpy(clone->triangles_sharing_sides,
		       msh3trg->triangles_sharing_sides,
		       3 * clone->N_triangles *
		       sizeof(*(clone->triangles_sharing_sides)));
	}
	clone->N_input_vertices = msh3trg->N_input_vertices;
	clone->input_vertices = malloc(clone->N_input_vertices *
				       sizeof(*(clone->input_vertices)));
	memcpy(clone->input_vertices,
	       msh3trg->input_vertices,
	       clone->N_input_vertices * sizeof(*(clone->input_vertices)));

	clone->N_input_segments = msh3trg->N_input_segments;
	if (clone->N_input_segments > 0) {
		clone->N_subsgm_x_inputsgm = 
			malloc(clone->N_input_segments *
			       sizeof(*(clone->N_subsgm_x_inputsgm)));
		memcpy(clone->N_subsgm_x_inputsgm,
		       msh3trg->N_subsgm_x_inputsgm,
		       clone->N_input_segments *
		       sizeof(*(clone->N_subsgm_x_inputsgm)));
		clone->meshvtx_x_inputsgm =
			calloc(clone->N_input_segments,
			       sizeof(*(clone->meshvtx_x_inputsgm)));
		for (uint32_t i = 0; i < clone->N_input_segments; i++) {
			if (clone->N_subsgm_x_inputsgm[i] == 0)
				continue;
			clone->meshvtx_x_inputsgm[i] =
				malloc((clone->N_subsgm_x_inputsgm[i]+1) *
				       sizeof(*(clone->meshvtx_x_inputsgm[i])));
			memcpy(clone->meshvtx_x_inputsgm[i],
			       msh3trg->meshvtx_x_inputsgm[i],
			       (clone->N_subsgm_x_inputsgm[i]+1) *
			       sizeof(*(clone->meshvtx_x_inputsgm[i])));
		}
	}    
	return clone;
}

void vcn_msh3trg_clear(vcn_msh3trg_t* msh3trg)
{
	if (msh3trg->N_vertices > 0) 
		free(msh3trg->vertices);
	if (msh3trg->N_edges > 0)
		free(msh3trg->edges);
	if (msh3trg->N_triangles > 0) {
		free(msh3trg->vertices_forming_triangles);
		if (NULL != msh3trg->triangles_sharing_sides)
			free(msh3trg->triangles_sharing_sides);
	}
	if (msh3trg->N_input_vertices > 0)
		free(msh3trg->input_vertices);
	if (msh3trg->N_input_segments > 0) {
		for (uint32_t i = 0; i < msh3trg->N_input_segments; i++) {
			if (msh3trg->N_subsgm_x_inputsgm[i] == 0)
				continue;
			free(msh3trg->meshvtx_x_inputsgm[i]);
		}
		free(msh3trg->N_subsgm_x_inputsgm);
		free(msh3trg->meshvtx_x_inputsgm);
	}
	memset(msh3trg, 0, sizeof(*msh3trg));
}

void vcn_msh3trg_destroy(vcn_msh3trg_t* msh3trg)
{
	vcn_msh3trg_finish(msh3trg);
	free(msh3trg);
}

bool nb_msh3trg_is_vtx_inside(const vcn_msh3trg_t *msh3trg,
			      const double x[2])
{
	bool is_inside = false;
	for (uint32_t i = 0; i < msh3trg->N_triangles; i++) {
		uint32_t id1 = msh3trg->vertices_forming_triangles[i * 3];
		uint32_t id2 = msh3trg->vertices_forming_triangles[i*3+1];
		uint32_t id3 = msh3trg->vertices_forming_triangles[i*3+2];
		double *v1 = &(msh3trg->vertices[id1*2]);
		double *v2 = &(msh3trg->vertices[id2*2]);
		double *v3 = &(msh3trg->vertices[id3*2]);
		if (vcn_utils2D_pnt_lies_in_trg(v1, v2, v3, x)) {
			is_inside = true;
			break;
		}
	}
	return is_inside;
}

vcn_graph_t* vcn_msh3trg_create_vtx_graph
                (const vcn_msh3trg_t *const restrict msh3trg)
{
	vcn_graph_t *graph = calloc(1, sizeof(*graph));
	graph->N = msh3trg->N_vertices;
	graph->N_adj = calloc(graph->N, sizeof(*(graph->N_adj)));
	graph->adj = calloc(graph->N, sizeof(*(graph->adj)));

	/* Connectivity Matrix stored in container */
	nb_container_t** l_nodal_CM = malloc(graph->N * sizeof(*l_nodal_CM));
	for (uint32_t i = 0; i < graph->N; i++)
		l_nodal_CM[i] = nb_container_create(NB_SORTED);
  
	for (uint32_t k = 0; k < msh3trg->N_triangles; k++) {
		for (uint32_t i = 0; i < 2; i++) {
			for (uint32_t j = i+1; j < 3; j++) {
				uint32_t* inode = malloc(sizeof(*inode));
				uint32_t* jnode = malloc(sizeof(*jnode));
				*inode = msh3trg->vertices_forming_triangles[k * 3 + i];
				*jnode = msh3trg->vertices_forming_triangles[k * 3 + j];
				uint32_t length1 = nb_container_get_length(l_nodal_CM[inode[0]]);
				uint32_t length2 = nb_container_get_length(l_nodal_CM[jnode[0]]);
				if (nb_container_exist(l_nodal_CM[inode[0]], jnode) == NULL)
					nb_container_insert(l_nodal_CM[inode[0]], jnode);
				if (nb_container_exist(l_nodal_CM[jnode[0]], inode) == NULL)
					nb_container_insert(l_nodal_CM[jnode[0]], inode);
				bool free_j = false;
				if (length1 == nb_container_get_length(l_nodal_CM[inode[0]]))
					free_j = true;
				if (length2 == nb_container_get_length(l_nodal_CM[jnode[0]]))
					free(inode);
				if (free_j) 
					free(jnode);
			}
		}
	}
	for (uint32_t i = 0; i < graph->N; i++) {
		graph->N_adj[i] = nb_container_get_length(l_nodal_CM[i]);
		graph->adj[i] = malloc(graph->N_adj[i] * sizeof(uint32_t));
    
		int j = 0;
		while (nb_container_is_not_empty(l_nodal_CM[i])) {
			uint32_t* inode = 
				nb_container_delete_first(l_nodal_CM[i]);
			graph->adj[i][j++] = *inode;
			free(inode);
		}
		nb_container_destroy(l_nodal_CM[i]);
	}
	free(l_nodal_CM);
	return graph;
}

vcn_graph_t* vcn_msh3trg_create_elem_graph
                (const vcn_msh3trg_t *const restrict msh3trg)
{
	if (NULL == msh3trg->triangles_sharing_sides)
		return NULL;

	vcn_graph_t* graph = calloc(1, sizeof(vcn_graph_t));

	graph->N = msh3trg->N_triangles;
	graph->N_adj = (uint32_t*) calloc(graph->N, sizeof(uint32_t));
	graph->adj = (uint32_t**) malloc(graph->N * sizeof(uint32_t*));
	for (uint32_t i = 0; i < graph->N; i++) {
		graph->N_adj[i] = 0;
		if (msh3trg->triangles_sharing_sides[i * 3] <
		    msh3trg->N_triangles)
			graph->N_adj[i] += 1;      

		if (msh3trg->triangles_sharing_sides[i*3+1] <
		    msh3trg->N_triangles)
			graph->N_adj[i] += 1;

		if (msh3trg->triangles_sharing_sides[i*3+2] <
		    msh3trg->N_triangles)
			graph->N_adj[i] += 1;
  }

	for (uint32_t i = 0; i < graph->N; i++) {
		graph->adj[i] = malloc(graph->N_adj[i] * sizeof(uint32_t));
		uint32_t cnt = 0;
		if (msh3trg->triangles_sharing_sides[i * 3] <
		    msh3trg->N_triangles)
			graph->adj[i][cnt++] = 
				msh3trg->triangles_sharing_sides[i * 3];

		if (msh3trg->triangles_sharing_sides[i*3+1] <
		    msh3trg->N_triangles)
			graph->adj[i][cnt++] =
				msh3trg->triangles_sharing_sides[i*3+1];
		
		if (msh3trg->triangles_sharing_sides[i*3+2] <
		    msh3trg->N_triangles)
			graph->adj[i][cnt++] =
				msh3trg->triangles_sharing_sides[i*3+2];
	}
	return graph;
}

void vcn_msh3trg_load_from_mesh(vcn_msh3trg_t *msh3trg,
				const vcn_mesh_t *const mesh)
{
	nb_trg_exporter_interface_t exp;
	exp.structure = msh3trg;
	exp.set_N_vtx = msh3trg_set_N_vtx;
	exp.malloc_vtx = msh3trg_malloc_vtx;
	exp.set_vtx = msh3trg_set_vtx;

	exp.set_N_edg = msh3trg_set_N_edg;
	exp.malloc_edg = msh3trg_malloc_edg;
	exp.set_edg = msh3trg_set_edg;

	exp.set_N_trg = msh3trg_set_N_trg;
	exp.malloc_trg = msh3trg_malloc_trg;
	exp.set_trg = msh3trg_set_trg;
	exp.set_trg_neighbours = msh3trg_set_trg_neighbours;

	exp.set_N_input_vtx = msh3trg_set_N_input_vtx;
	exp.malloc_input_vtx = msh3trg_malloc_input_vtx;
	exp.set_input_vtx = msh3trg_set_input_vtx;

	exp.set_N_input_sgm = msh3trg_set_N_input_sgm;
	exp.malloc_input_sgm_table = msh3trg_malloc_input_sgm_table;
	exp.input_sgm_set_N_vtx = msh3trg_input_sgm_set_N_vtx;
	exp.input_sgm_malloc_vtx = msh3trg_input_sgm_malloc_vtx;
	exp.input_sgm_set_vtx = msh3trg_input_sgm_set_vtx;

	vcn_mesh_export(mesh, &exp);
}

static void msh3trg_set_N_vtx(void *msh3trg_ptr, uint32_t N)
{
	vcn_msh3trg_t *msh3trg = msh3trg_ptr;
	msh3trg->N_vertices = N;
}

static void msh3trg_malloc_vtx(void *msh3trg_ptr)
{
	vcn_msh3trg_t *msh3trg = msh3trg_ptr;
	uint32_t memsize = 2 * msh3trg->N_vertices * 
		sizeof(*(msh3trg->vertices));
	msh3trg->vertices = malloc(memsize);
}

static void msh3trg_set_vtx(void *msh3trg_ptr, uint32_t i,
			    double x, double y)
{
	vcn_msh3trg_t *msh3trg = msh3trg_ptr;
	msh3trg->vertices[i * 2] = x;
	msh3trg->vertices[i*2+1] = y;
}

static void msh3trg_set_N_edg(void *msh3trg_ptr, uint32_t N)
{
	vcn_msh3trg_t *msh3trg = msh3trg_ptr;
	msh3trg->N_edges = N;
}

static void msh3trg_malloc_edg(void *msh3trg_ptr)
{
	vcn_msh3trg_t *msh3trg = msh3trg_ptr;
	uint32_t memsize = 2 * msh3trg->N_edges * sizeof(*(msh3trg->edges));
	msh3trg->edges = malloc(memsize);
}

static void msh3trg_set_edg(void *msh3trg_ptr, uint32_t i,
			    uint32_t v1, uint32_t v2)
{
	vcn_msh3trg_t *msh3trg = msh3trg_ptr;
	msh3trg->edges[i * 2] = v1;
	msh3trg->edges[i*2+1] = v2;
}

static void msh3trg_set_N_trg(void *msh3trg_ptr, uint32_t N)
{
	vcn_msh3trg_t *msh3trg = msh3trg_ptr;
	msh3trg->N_triangles = N;
}

static void msh3trg_malloc_trg(void *msh3trg_ptr, bool include_neighbours)
{
	vcn_msh3trg_t *msh3trg = msh3trg_ptr;
	uint32_t vtx_memsize = 3 * msh3trg->N_triangles *
		sizeof(*(msh3trg->vertices_forming_triangles));
	msh3trg->vertices_forming_triangles = malloc(vtx_memsize);

	if (include_neighbours) {
		uint32_t nbg_memsize = 3 * msh3trg->N_triangles *
			sizeof(*(msh3trg->triangles_sharing_sides));
		msh3trg->triangles_sharing_sides = malloc(nbg_memsize);
	}
}

static void msh3trg_set_trg(void *msh3trg_ptr, uint32_t i,
			    uint32_t v1, uint32_t v2, uint32_t v3)
{
	vcn_msh3trg_t *msh3trg = msh3trg_ptr;
	msh3trg->vertices_forming_triangles[i * 3] = v1;
	msh3trg->vertices_forming_triangles[i*3+1] = v2;
	msh3trg->vertices_forming_triangles[i*3+2] = v3;
}

static void msh3trg_set_trg_neighbours(void *msh3trg_ptr,
				       uint32_t i, uint32_t t1,
				       uint32_t t2, uint32_t t3)

{
	vcn_msh3trg_t *msh3trg = msh3trg_ptr;
	msh3trg->triangles_sharing_sides[i * 3] = t1;
	msh3trg->triangles_sharing_sides[i*3+1] = t2;
	msh3trg->triangles_sharing_sides[i*3+2] = t3;
}

static void msh3trg_set_N_input_vtx(void *msh3trg_ptr, uint32_t N)
{
	vcn_msh3trg_t *msh3trg = msh3trg_ptr;
	msh3trg->N_input_vertices = N;
}

static void msh3trg_malloc_input_vtx(void *msh3trg_ptr)
{
	vcn_msh3trg_t *msh3trg = msh3trg_ptr;
	uint32_t memsize = msh3trg->N_input_vertices *
		sizeof(*(msh3trg->input_vertices));
	msh3trg->input_vertices = malloc(memsize);
}

static void msh3trg_set_input_vtx(void *msh3trg_ptr, uint32_t i,
				  uint32_t vtx_id)
{
	vcn_msh3trg_t *msh3trg = msh3trg_ptr;
	msh3trg->input_vertices[i] = vtx_id;
}

static void msh3trg_set_N_input_sgm(void *msh3trg_ptr, uint32_t N)
{
	vcn_msh3trg_t *msh3trg = msh3trg_ptr;
	msh3trg->N_input_segments = N;
}

static void msh3trg_malloc_input_sgm_table(void *msh3trg_ptr)
{
	vcn_msh3trg_t *msh3trg = msh3trg_ptr;
	uint32_t sgm_memsize = msh3trg->N_input_segments *
		sizeof(*(msh3trg->N_subsgm_x_inputsgm));
	msh3trg->N_subsgm_x_inputsgm = malloc(sgm_memsize);
	uint32_t vtx_memsize = msh3trg->N_input_segments *
		sizeof(*(msh3trg->meshvtx_x_inputsgm));
	msh3trg->meshvtx_x_inputsgm = malloc(vtx_memsize);
}

static void msh3trg_input_sgm_set_N_vtx(void *msh3trg_ptr, uint32_t i,
					uint32_t N)
{
	vcn_msh3trg_t *msh3trg = msh3trg_ptr;
	msh3trg->N_subsgm_x_inputsgm[i] = N;
}

static void msh3trg_input_sgm_malloc_vtx(void *msh3trg_ptr, uint32_t i)
{
	vcn_msh3trg_t *msh3trg = msh3trg_ptr;
	uint32_t memsize = (1 + msh3trg->N_subsgm_x_inputsgm[i]) *
		sizeof(**(msh3trg->meshvtx_x_inputsgm));
	msh3trg->meshvtx_x_inputsgm[i] = malloc(memsize);
}

static void msh3trg_input_sgm_set_vtx(void *msh3trg_ptr, uint32_t isgm,
				      uint32_t ivtx, uint32_t vtx_id)
{
	vcn_msh3trg_t *msh3trg = msh3trg_ptr;
	msh3trg->meshvtx_x_inputsgm[isgm][ivtx] = vtx_id;
}

void vcn_msh3trg_relabel(vcn_msh3trg_t* msh3trg,
			    uint32_t* (*labeling)(const vcn_graph_t *const))
{
	;/* PENDING */
}

void vcn_msh3trg_disable_single_point_connections
                      (const vcn_msh3trg_t *const restrict msh3trg,
		       bool* enabled_elements)
{
	/* Allocate lists to store triangles per vertex */
	uint32_t *N_trg_x_vtx = calloc(msh3trg->N_vertices, 
				       sizeof(*N_trg_x_vtx));
	uint32_t **trg_x_vtx = malloc(msh3trg->N_vertices *
				      sizeof(*trg_x_vtx));
	for (uint32_t i = 0; i < msh3trg->N_vertices; i++)
		trg_x_vtx[i] = calloc(10, sizeof(*(trg_x_vtx[i])));
      
	/* Iterate over triangles to found relations */
	for (uint32_t i = 0; i < msh3trg->N_triangles; i++) {
		if (!enabled_elements[i])
			continue;
		uint32_t v1 = msh3trg->vertices_forming_triangles[i * 3];
		uint32_t v2 = msh3trg->vertices_forming_triangles[i*3+1];
		uint32_t v3 = msh3trg->vertices_forming_triangles[i*3+2];
		trg_x_vtx[v1][N_trg_x_vtx[v1]++] = i;
		trg_x_vtx[v2][N_trg_x_vtx[v2]++] = i;
		trg_x_vtx[v3][N_trg_x_vtx[v3]++] = i;
	}

	/* Detect one point connections */
	for (uint32_t i = 0; i < msh3trg->N_vertices; i++) {
		if(N_trg_x_vtx[i] < 2)
			continue;

		uint32_t itrg = trg_x_vtx[i][0];

		uint32_t itrg_twist = itrg;
		bool twist_around = false;

		while (itrg_get_right_triangle(msh3trg, 
					       enabled_elements,
					       itrg_twist, i) !=  msh3trg->N_triangles) {
			itrg_twist = itrg_get_right_triangle(msh3trg, enabled_elements,
							     itrg_twist, i);
			if (itrg_twist == itrg) {
				twist_around = true;
				break;
			}
		}
		if (itrg_twist == itrg && twist_around)
			continue;

		uint32_t N_trg_fan = 0;
		uint32_t trg_fan[12];
		do {
			trg_fan[N_trg_fan++] = itrg_twist;
			itrg_twist = itrg_get_left_triangle(msh3trg, enabled_elements,
							    itrg_twist, i);
		} while (itrg_twist != msh3trg->N_triangles);

		if (N_trg_fan == N_trg_x_vtx[i])
			continue;

		for (uint32_t j = 0; j < N_trg_x_vtx[i]; j++) {
			bool is_inside_the_fan = false;
			for (uint32_t k = 0; k < N_trg_fan; k++) {
				if (trg_x_vtx[i][j] == trg_fan[k]) {
					is_inside_the_fan = true;
					break;
				}
			}
			if (is_inside_the_fan)
				continue;
			enabled_elements[trg_x_vtx[i][j]] = false;
		}
	}
  
	/* Free memory */
	for (uint32_t i = 0; i < msh3trg->N_vertices; i++)
		free(trg_x_vtx[i]);
	free(trg_x_vtx);
	free(N_trg_x_vtx);
}

static inline uint32_t itrg_get_right_triangle
                         (const vcn_msh3trg_t *const delaunay, 
			  const bool *const enabled_elements,
			  uint32_t itrg, uint32_t ivtx)
{
	if (delaunay->vertices_forming_triangles[itrg * 3] == ivtx) {
		uint32_t id = delaunay->triangles_sharing_sides[itrg * 3];
		if (id < delaunay->N_triangles)
			if (enabled_elements[id])
				return id;
		return delaunay->N_triangles;
	}
	if (delaunay->vertices_forming_triangles[itrg*3+1] == ivtx) {
		uint32_t id = delaunay->triangles_sharing_sides[itrg*3+1];
		if (id < delaunay->N_triangles)
			if (enabled_elements[id])
				return id;
		return delaunay->N_triangles;
	}
	if (delaunay->vertices_forming_triangles[itrg*3+2] == ivtx) {
		uint32_t id = delaunay->triangles_sharing_sides[itrg*3+2];
		if (id < delaunay->N_triangles)
			if (enabled_elements[id])
				return id;
		return delaunay->N_triangles;
	}
	return delaunay->N_triangles;
}

static inline uint32_t itrg_get_left_triangle
                         (const vcn_msh3trg_t *const delaunay,
			  const bool *const enabled_elements,
			  uint32_t itrg, uint32_t ivtx)
{
	if(delaunay->vertices_forming_triangles[itrg * 3] == ivtx) {
		uint32_t id = delaunay->triangles_sharing_sides[itrg*3+2];
		if(id < delaunay->N_triangles)
			if(enabled_elements[id])
				return id;
		return delaunay->N_triangles;
	}
	if (delaunay->vertices_forming_triangles[itrg*3+1] == ivtx) {
		uint32_t id = delaunay->triangles_sharing_sides[itrg * 3];
		if(id < delaunay->N_triangles)
			if(enabled_elements[id])
				return id;
		return delaunay->N_triangles;
	}
	if (delaunay->vertices_forming_triangles[itrg*3+2] == ivtx) {
		uint32_t id = delaunay->triangles_sharing_sides[itrg*3+1];
		if(id < delaunay->N_triangles)
			if(enabled_elements[id])
				return id;
		return delaunay->N_triangles;
	}
	return delaunay->N_triangles;
}
