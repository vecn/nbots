#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>

#include "vcn/container_bot/container.h"
#include "vcn/container_bot/iterator.h"
#include "vcn/geometric_bot/knn/bins2D.h"
#include "vcn/geometric_bot/knn/bins2D_iterator.h"
#include "vcn/geometric_bot/mesh/elements2D/triangles_struct.h"
#include "vcn/geometric_bot/mesh/elements2D/triangles.h"

#include "../mesh2D_structs.h"

static uint32_t key_uint(const void *const uint_ptr);
static uint32_t itrg_get_right_triangle
                         (const vcn_msh3trg_t *const delaunay, 
			  const bool *const enabled_elements,
			  uint32_t itrg, uint32_t ivtx);
static uint32_t itrg_get_left_triangle
                         (const vcn_msh3trg_t *const delaunay,
			  const bool *const enabled_elements,
			  uint32_t itrg, uint32_t ivtx);
static void mesh_2_msh3trg_cast_vertices(const vcn_mesh_t *const mesh,
					 vcn_msh3trg_t *msh3trg);
static void mesh_2_msh3trg_cast_edges(const vcn_mesh_t *const mesh,
				      vcn_msh3trg_t *msh3trg);
static void mesh_2_msh3trg_cast_trg_neighbours(const vcn_mesh_t *const mesh,
					       vcn_msh3trg_t *msh3trg);
static void mesh_2_msh3trg_cast_trg_and_alloc_id(vcn_mesh_t *mesh,
						 vcn_msh3trg_t *msh3trg,
						 bool include_neighbours);
static void mesh_2_msh3trg_cast_input_vtx(const vcn_mesh_t *const mesh,
					  vcn_msh3trg_t *msh3trg);
static void mesh_2_msh3trg_cast_input_sgm(const vcn_mesh_t *const mesh,
					  vcn_msh3trg_t * msh3trg);
static void mesh_2_msh3trg_free_trg_ids(vcn_mesh_t *mesh);
static void mesh_alloc_vtx_ids(vcn_mesh_t *mesh);
static void mesh_free_vtx_ids(vcn_mesh_t *mesh);

inline vcn_msh3trg_t* vcn_msh3trg_create(void)
{
	vcn_msh3trg_t* delaunay = calloc(1, sizeof(*delaunay));
	return delaunay;
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
	vcn_msh3trg_clear(msh3trg);
	free(msh3trg);
}

vcn_graph_t* vcn_msh3trg_create_vtx_graph
                (const vcn_msh3trg_t *const restrict msh3trg)
{
	vcn_graph_t *graph = calloc(1, sizeof(*graph));
	graph->N = msh3trg->N_vertices;
	graph->N_adj = calloc(graph->N, sizeof(uint32_t));
	graph->adj = calloc(graph->N, sizeof(uint32_t*));

	/* Connectivity Matrix stored in container */
	vcn_container_t** l_nodal_CM = malloc(graph->N * sizeof(*l_nodal_CM));
	for (uint32_t i = 0; i < graph->N; i++) {
		l_nodal_CM[i] = vcn_container_create(VCN_CONTAINER_SORTED);
		vcn_container_set_key_generator(l_nodal_CM[i], key_uint);
	}
  
	for (uint32_t k = 0; k < msh3trg->N_triangles; k++) {
		for (uint32_t i = 0; i < 2; i++) {
			for (uint32_t j = i+1; j < 3; j++) {
				uint32_t* inode = malloc(sizeof(*inode));
				uint32_t* jnode = malloc(sizeof(*jnode));
				*inode = msh3trg->vertices_forming_triangles[k * 3 + i];
				*jnode = msh3trg->vertices_forming_triangles[k * 3 + j];
				uint32_t length1 = vcn_container_get_length(l_nodal_CM[inode[0]]);
				uint32_t length2 = vcn_container_get_length(l_nodal_CM[jnode[0]]);
				if (vcn_container_exist(l_nodal_CM[inode[0]], jnode) == NULL)
					vcn_container_insert(l_nodal_CM[inode[0]], jnode);
				if (vcn_container_exist(l_nodal_CM[jnode[0]], inode) == NULL)
					vcn_container_insert(l_nodal_CM[jnode[0]], inode);
				bool free_j = false;
				if (length1 == vcn_container_get_length(l_nodal_CM[inode[0]]))
					free_j = true;
				if (length2 == vcn_container_get_length(l_nodal_CM[jnode[0]]))
					free(inode);
				if (free_j) 
					free(jnode);
			}
		}
	}
	for (uint32_t i = 0; i < graph->N; i++) {
		graph->N_adj[i] = vcn_container_get_length(l_nodal_CM[i]);
		graph->adj[i] = malloc(graph->N_adj[i] * sizeof(uint32_t));
    
		int j = 0;
		while (vcn_container_is_not_empty(l_nodal_CM[i])) {
			uint32_t* inode = 
				vcn_container_delete_first(l_nodal_CM[i]);
			graph->adj[i][j++] = *inode;
			free(inode);
		}
		vcn_container_destroy(l_nodal_CM[i]);
	}
	free(l_nodal_CM);
	return graph;
}

static inline uint32_t key_uint(const void *const uint_ptr)
{
	return *((uint32_t*)uint_ptr);
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

vcn_msh3trg_t* vcn_mesh_get_msh3trg
                 (const vcn_mesh_t *const mesh,
		  bool include_edges,
		  bool include_triangles,
		  bool include_neighbours,
		  bool include_input_vertices,
		  bool include_input_segments)
{
	vcn_msh3trg_t* msh3trg = vcn_msh3trg_create();
	mesh_alloc_vtx_ids((vcn_mesh_t*)mesh);

	mesh_2_msh3trg_cast_vertices(mesh, msh3trg);

	if (include_edges)
		mesh_2_msh3trg_cast_edges(mesh, msh3trg);

	if(include_triangles && vcn_container_is_not_empty(mesh->ht_trg))
		mesh_2_msh3trg_cast_trg_and_alloc_id((vcn_mesh_t*)mesh, msh3trg, 
						     include_neighbours);


	if (include_input_vertices)
		mesh_2_msh3trg_cast_input_vtx(mesh, msh3trg);

	if (include_input_segments && mesh->N_input_sgm > 0)
		mesh_2_msh3trg_cast_input_sgm(mesh, msh3trg);

	/* Free memory */
	mesh_free_vtx_ids((vcn_mesh_t*)mesh); /* Casting to non-const */

	if (include_neighbours)
		mesh_2_msh3trg_free_trg_ids((vcn_mesh_t*)mesh); /* Casting to non-const */

	/* Return data */
	return msh3trg;
}

void vcn_mesh_trg2D_relabel(vcn_msh3trg_t* msh3trg,
			    uint32_t* (*labeling)(const vcn_graph_t *const))
{

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


static void mesh_2_msh3trg_cast_vertices(const vcn_mesh_t *const restrict mesh,
					 vcn_msh3trg_t * restrict msh3trg)
{
	msh3trg->N_vertices = vcn_bins2D_get_length(mesh->ug_vtx);
	msh3trg->vertices =
		malloc(2 * msh3trg->N_vertices * sizeof(*(msh3trg->vertices)));
	vcn_bins2D_iter_t* iter = vcn_bins2D_iter_create();
	vcn_bins2D_iter_set_bins(iter, mesh->ug_vtx);
	while (vcn_bins2D_iter_has_more(iter)) {
		const msh_vtx_t* vtx = vcn_bins2D_iter_get_next(iter);
		uint32_t id = *(uint32_t*)((void**)vtx->attr)[0];
		msh3trg->vertices[id * 2] = vtx->x[0]/mesh->scale + mesh->xdisp;
		msh3trg->vertices[id*2+1] = vtx->x[1]/mesh->scale + mesh->ydisp;
	}
	vcn_bins2D_iter_destroy(iter);
}

static void mesh_2_msh3trg_cast_edges(const vcn_mesh_t *const restrict mesh,
				      vcn_msh3trg_t * restrict msh3trg)
{
	msh3trg->N_edges = vcn_container_get_length(mesh->ht_edge);
	msh3trg->edges = malloc(2 * msh3trg->N_edges * sizeof(*(msh3trg->edges)));

	uint32_t i = 0;
	vcn_iterator_t *iter = vcn_iterator_create();
	vcn_iterator_set_container(iter, mesh->ht_edge);
	while (vcn_iterator_has_more(iter)) {
		msh_edge_t *edge = (msh_edge_t*) vcn_iterator_get_next(iter);
		msh3trg->edges[i * 2] = *(uint32_t*)((void**)edge->v1->attr)[0];
		msh3trg->edges[i*2+1] = *(uint32_t*)((void**)edge->v2->attr)[0];
		i += 1;
	}
	vcn_iterator_destroy(iter);
}

static void mesh_2_msh3trg_cast_trg_neighbours(const vcn_mesh_t *const restrict mesh,
					       vcn_msh3trg_t * restrict msh3trg)
{
	msh3trg->triangles_sharing_sides = 
		malloc(3 * msh3trg->N_triangles *
		       sizeof(*(msh3trg->triangles_sharing_sides)));

    vcn_iterator_t* trg_iter =  vcn_iterator_create();
    vcn_iterator_set_container(trg_iter, mesh->ht_trg);
    while (vcn_iterator_has_more(trg_iter)) {
	    msh_trg_t* trg = (msh_trg_t*)vcn_iterator_get_next(trg_iter);
	    uint32_t id = ((uint32_t*)((void**)trg->attr)[0])[0];
	    uint32_t id1 = msh3trg->N_triangles;
	    if (NULL != trg->t1)
		    id1 = ((uint32_t*)((void**)trg->t1->attr)[0])[0];
	    uint32_t id2 = msh3trg->N_triangles;
	    if (NULL != trg->t2)
		    id2 = ((uint32_t*)((void**)trg->t2->attr)[0])[0];
	    uint32_t id3 = msh3trg->N_triangles;
	    if (NULL != trg->t3)
		    id3 = ((uint32_t*)((void**)trg->t3->attr)[0])[0];

	    msh3trg->triangles_sharing_sides[id * 3] = id1;
	    msh3trg->triangles_sharing_sides[id*3+1] = id2;
	    msh3trg->triangles_sharing_sides[id*3+2] = id3;
    }
    vcn_iterator_destroy(trg_iter);
}

static void mesh_2_msh3trg_cast_trg_and_alloc_id
                                   (vcn_mesh_t * restrict mesh,
				    vcn_msh3trg_t * restrict msh3trg,
				    bool include_neighbours)
{
	msh3trg->N_triangles = vcn_container_get_length(mesh->ht_trg);
	msh3trg->vertices_forming_triangles =
		malloc(3 * msh3trg->N_triangles * 
		       sizeof(*(msh3trg->vertices_forming_triangles)));

	vcn_iterator_t* trg_iter = vcn_iterator_create();
	vcn_iterator_set_container(trg_iter, mesh->ht_trg);
	uint32_t i = 0;
	while (vcn_iterator_has_more(trg_iter)) {
		msh_trg_t* trg = (msh_trg_t*)vcn_iterator_get_next(trg_iter);
		msh3trg->vertices_forming_triangles[i * 3] =
			*((uint32_t*)((void**)trg->v1->attr)[0]);
		msh3trg->vertices_forming_triangles[i*3+1] = 
			*((uint32_t*)((void**)trg->v2->attr)[0]);
		msh3trg->vertices_forming_triangles[i*3+2] =
			*((uint32_t*)((void**)trg->v3->attr)[0]);
		if (include_neighbours) {
			void** attr = malloc(2 * sizeof(*attr));
			uint32_t* id = malloc(sizeof(*id));
			id[0] = i;
			attr[0] = id;
			attr[1] = trg->attr;
			trg->attr = attr;
		}
		i++;
	}
	vcn_iterator_destroy(trg_iter);
  
	if(include_neighbours)
		mesh_2_msh3trg_cast_trg_neighbours(mesh, msh3trg);
}

static void mesh_2_msh3trg_cast_input_vtx(const vcn_mesh_t *const restrict mesh,
					  vcn_msh3trg_t * restrict msh3trg)
{
	msh3trg->N_input_vertices = mesh->N_input_vtx;
	msh3trg->input_vertices = 
		malloc(msh3trg->N_input_vertices * sizeof(*(msh3trg->input_vertices)));
	for (uint32_t i = 0; i < msh3trg->N_input_vertices; i++) {
		if (NULL == mesh->input_vtx[i]) {
			msh3trg->input_vertices[i] = msh3trg->N_vertices;
			continue;
		}
		void** attr = mesh->input_vtx[i]->attr;
		uint32_t* id = (uint32_t*)attr[0];
		msh3trg->input_vertices[i] = id[0];
	}
}

static void mesh_2_msh3trg_cast_input_sgm(const vcn_mesh_t *const restrict mesh,
					  vcn_msh3trg_t * restrict msh3trg)
{
	msh3trg->N_input_segments = mesh->N_input_sgm;
	msh3trg->N_subsgm_x_inputsgm =
		calloc(msh3trg->N_input_segments, sizeof(*(msh3trg->N_subsgm_x_inputsgm)));
	msh3trg->meshvtx_x_inputsgm = 
		calloc(msh3trg->N_input_segments,  sizeof(*(msh3trg->meshvtx_x_inputsgm)));

	for (uint32_t i = 0; i < msh3trg->N_input_segments; i++) {
		msh_edge_t* sgm = mesh->input_sgm[i];
		uint32_t counter = 0;
		while (NULL != sgm) {
			counter++;
			sgm = medge_subsgm_next(sgm);
		}
		msh3trg->N_subsgm_x_inputsgm[i] = counter;
		if (counter > 0)
			msh3trg->meshvtx_x_inputsgm[i] =
				calloc(counter+1, sizeof(*(msh3trg->meshvtx_x_inputsgm[i])));
  }

	for (uint32_t i=0; i < msh3trg->N_input_segments; i++) {
		if (NULL == mesh->input_sgm[i])
			continue;
		msh_edge_t* sgm = mesh->input_sgm[i];
		msh_edge_t* sgm_prev = sgm;
		sgm = medge_subsgm_next(sgm);
		if (NULL == sgm) {
			msh3trg->meshvtx_x_inputsgm[i][0] = 
				((uint32_t*)((void**)sgm_prev->v1->attr)[0])[0];
			msh3trg->meshvtx_x_inputsgm[i][1] =
				((uint32_t*)((void**)sgm_prev->v2->attr)[0])[0];
		} else {
			uint32_t idx = 0;
			uint32_t id_chain;
			uint32_t id1 = ((uint32_t*)((void**)sgm_prev->v1->attr)[0])[0];
			uint32_t id2 = ((uint32_t*)((void**)sgm_prev->v2->attr)[0])[0];
			uint32_t id1n = ((uint32_t*)((void**)sgm->v1->attr)[0])[0];
			uint32_t id2n = ((uint32_t*)((void**)sgm->v2->attr)[0])[0];
			if (id2 == id1n || id2 == id2n) {
				msh3trg->meshvtx_x_inputsgm[i][idx++] =  id1;
				msh3trg->meshvtx_x_inputsgm[i][idx++] =  id2;
				id_chain = id2;
			} else {
				msh3trg->meshvtx_x_inputsgm[i][idx++] =  id2;
				msh3trg->meshvtx_x_inputsgm[i][idx++] =  id1;
				id_chain = id1;
			}
			while (NULL != sgm) {
				sgm_prev = sgm;
				uint32_t id1 = ((uint32_t*)((void**)sgm_prev->v1->attr)[0])[0];
				uint32_t id2 = ((uint32_t*)((void**)sgm_prev->v2->attr)[0])[0];
				if (id1 == id_chain) {
					msh3trg->meshvtx_x_inputsgm[i][idx++] =  id2;
					id_chain = id2;
				} else {
					msh3trg->meshvtx_x_inputsgm[i][idx++] =  id1;
					id_chain = id1;
				}
				sgm = medge_subsgm_next(sgm);
			}
		}
	}
}

static void mesh_2_msh3trg_free_trg_ids(vcn_mesh_t *restrict mesh)
{
	vcn_iterator_t* iter = vcn_iterator_create();
	vcn_iterator_set_container(iter, mesh->ht_trg);
	while (vcn_iterator_has_more(iter)) {
		msh_trg_t* trg = (msh_trg_t*)vcn_iterator_get_next(iter);
		void** attr = (void**)trg->attr;
		trg->attr = attr[1];
		free(attr[0]);
		free(attr);
	}
	vcn_iterator_destroy(iter);
}

static void mesh_alloc_vtx_ids(vcn_mesh_t * restrict mesh)
{
	vcn_bins2D_iter_t* iter = vcn_bins2D_iter_create();
	vcn_bins2D_iter_set_bins(iter, mesh->ug_vtx);
	int i = 0;
	while (vcn_bins2D_iter_has_more(iter)) {
		msh_vtx_t* vtx = (msh_vtx_t*) vcn_bins2D_iter_get_next(iter);
		void** attr = malloc(2 * sizeof(*attr));
		uint32_t* id = malloc(sizeof(*id));
		id[0] = i++;
		attr[0] = id;
		attr[1] = vtx->attr;
		vtx->attr = attr;
	}
	vcn_bins2D_iter_destroy(iter);
}

static void mesh_free_vtx_ids(vcn_mesh_t *mesh)
{
	vcn_bins2D_iter_t* iter = vcn_bins2D_iter_create();
	vcn_bins2D_iter_set_bins(iter, mesh->ug_vtx);
	while (vcn_bins2D_iter_has_more(iter)) {
		msh_vtx_t* vtx = (msh_vtx_t*) vcn_bins2D_iter_get_next(iter);
		void** attr = vtx->attr;
		vtx->attr = attr[1];
		free(attr[0]);
		free(attr);
	}
	vcn_bins2D_iter_destroy(iter);
}
