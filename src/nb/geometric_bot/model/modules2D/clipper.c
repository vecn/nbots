#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>
#include <math.h>

#include "nb/memory_bot.h"
#include "nb/container_bot/container.h"
#include "nb/container_bot/iterator.h"
#include "nb/geometric_bot/utils2D.h"
#include "nb/geometric_bot/mesh/mesh2D.h"
#include "nb/geometric_bot/mesh/partition.h"
#include "nb/geometric_bot/mesh/modules2D/area_analizer.h"

#include "vtx.h"
#include "edge.h"
#include "nb/geometric_bot/model/model2D_struct.h"
#include "nb/geometric_bot/model/model2D.h"
#include "nb/geometric_bot/model/modules2D/verifier.h"
#include "nb/geometric_bot/model/modules2D/clipper.h"

#define MIN_LENGTH_X_SGM 1e-4

#define GET_PVTX(model, i) (&((model)->vertex[(i)*2]))
#define GET_1_EDGE_VTX(model, i) ((model)->edge[(i) * 2])
#define GET_2_EDGE_VTX(model, i) ((model)->edge[(i)*2+1])
#define MAX(a, b) (((a)>(b))?(a):(b))

typedef struct {
	double scale;
	double xdisp;
	double ydisp;
} scaling_data;

typedef struct {
	nb_intersect_t status;
	edge_t *sgm1;
	edge_t *sgm2;
	double intersection[2];
} ipack_t;

typedef struct {
	double length1;
	double length2;
	double dist_1A_2A;
	double dist_1A_2B;
	double dist_1B_2A;
	double dist_1B_2B;
	double max_dist;
} collinear_dists;

static void get_combined_elements(const vcn_model_t *model1,
				  const vcn_model_t *model2,
				  nb_container_t *ht_vtx,
				  nb_container_t *avl_sgm);
static void search_intersections_in_edges(const nb_container_t *edges,
					  nb_container_t *intersections);
static void set_combined_model(vcn_model_t *model,
			       nb_container_t *ht_vtx,
			       nb_container_t *avl_sgm,
			       const vcn_model_t *model1,
			       const vcn_model_t *model2);
static void set_vertices(vcn_model_t *model,
			 nb_container_t *ht_vtx);
static void set_segments(vcn_model_t *model,
			 nb_container_t *avl_sgm);
static void set_holes(const vcn_model_t *model1,
		      const vcn_model_t *model2,
		      vcn_model_t *model);
static uint32_t mask_true_centroids(const vcn_model_t *model1,
				    const vcn_model_t *model2,
				    char *mask, const double *centroids,
				    uint32_t N_centroids);
static void rescale_model(vcn_model_t *model, const scaling_data *scaling);
static void* ipack_create(void);
static void ipack_destroy(void *ipack_ptr);
static int8_t ipack_comparer(const void *ipack1_ptr,
			     const void *ipack2_ptr);

static void get_scaling_data(const vcn_model_t *const model1,
			     const vcn_model_t *const model2,
			     scaling_data *scaling);
static void scale_and_displace(vcn_model_t *model, 
			       const scaling_data *scaling);
static void search_intersections(nb_container_t *intersections,
				 const vcn_model_t *model1,
				 const vcn_model_t *model2,
				 nb_container_t *vertices,
				 nb_container_t *edges);
static void insert_edges_and_vtx(const vcn_model_t *model,
				 nb_container_t *vertices,
				 nb_container_t *edges,
				 edge_t **edges_array,
				 bool check_if_edges_exist);
static vtx_t* insert_vtx_if_not_exist(const vcn_model_t *const model,
				      nb_container_t *vertices,
				      uint32_t vtx_id);
static edge_t* exist_edge_and_insert_if_not(nb_container_t* edges,
					    edge_t* edge);
static ipack_t* get_intersection_pack(const edge_t *const sgm1,
				      const edge_t *const sgm2);
static bool is_true_intersection(const edge_t *sgm1, const edge_t *sgm2,
				 nb_intersect_t status);
static bool if_collinear_sgm_are_not_intersected(const edge_t *sgm1,
						 const edge_t *sgm2);
static void process_segment_intersections(nb_container_t* avl_sgm,
					  nb_container_t* sgm_intersect,
					  nb_container_t* ht_vtx);
static void process_ipack(ipack_t* ipack,
			  nb_container_t* avl_sgm,
			  nb_container_t* sgm_intersect,
			  nb_container_t* ht_vtx);
static void process_intersected_sgm(ipack_t *ipack,
				    nb_container_t* avl_sgm,
				    nb_container_t* sgm_intersect,
				    nb_container_t* ht_vtx);
static void process_collinear_sgm(ipack_t *ipack,
				  nb_container_t* avl_sgm,
				  nb_container_t* sgm_intersect);
static void calculate_collinear_dists(const edge_t *sgm1,
				      const edge_t *sgm2,
				      collinear_dists *cdists);
static bool collinear_are_contained(collinear_dists *cdists);
static bool collinear_are_intersected(collinear_dists *cdists);
static void process_collinear_intersected(ipack_t *ipack,
					  collinear_dists *cdists,
					  nb_container_t* avl_sgm,
					  nb_container_t* sgm_intersect);
static void collinear_merge_into_sgm1(edge_t *sgm1, const edge_t *sgm2, 
				      const collinear_dists *cdists);
static void update_ipacks_sgm1(ipack_t *ipack,
			       nb_container_t *sgm_intersect,
			       nb_container_t *sgm_intersect_aux,
			       const edge_t *existing_sgm1);
static edge_t *ipack_get_opp_sgm(const ipack_t *ipack,
				 const edge_t *sgm);
static void update_ipacks_sgm2(ipack_t *ipack,
			       nb_container_t *sgm_intersect,
			       nb_container_t *sgm_intersect_aux,
			       const edge_t *existing_sgm1);
static void move_out_ipacks_from_aux(nb_container_t *sgm_intersect,
				     nb_container_t *sgm_intersect_aux);
static void process_collinear_contained(ipack_t *ipack,
					collinear_dists *cdists,
					nb_container_t* avl_sgm,
					nb_container_t* sgm_intersect);
static edge_t *get_smallest_sgm(ipack_t *ipack, collinear_dists *cdists);
static void remove_short_segments(nb_container_t* avl_sgm,
				  nb_container_t* ht_vtx,
				  double min_length_x_sgm);
static void get_short_sgm(const nb_container_t *avl_sgm,
			  nb_container_t *short_sgm,
			  double min_length_x_sgm);
static void collapse_sgm_into_v1(nb_container_t *avl_sgm,
				 nb_container_t *ht_vtx, edge_t *sgm);
static void get_connected_segments(const nb_container_t *avl_sgm,
				   const edge_t *sgm,
				   nb_container_t *list_sgm);
static void link_connected_segments_to_v1(const edge_t *sgm,
					  nb_container_t *avl_sgm,
					  nb_container_t *list_sgm);
static void v1_to_centroid_if_connections(const nb_container_t *list_sgm,
					  nb_container_t *ht_vtx,
					  edge_t *sgm);
static void update_segments_connected(nb_container_t *avl_sgm,
				      nb_container_t *list_sgm);
static void set_as_initial_vtx_in_edges(nb_container_t *edges);
static void delete_unused_vertices(nb_container_t *vertices);
static ipack_t* search_ipack_with_sgm
				(const nb_container_t *const packs,
				 const edge_t *const edge);
static ipack_t* search_ipack_with_both_sgm(const nb_container_t *const packs,
					   const edge_t *const edge1,
					   const edge_t *const edge2);
static void split_segment_by_vertex(nb_container_t* avl_sgm,
				    nb_container_t* sgm_intersect,
				    edge_t* sgm, vtx_t* vtx);
static void set_new_segments_ipacks(const edge_t *sgm,
				    const edge_t *new_sgm,
				    nb_container_t *sgm_intersect,
				    nb_container_t *sgm_intersect_aux,
				    const edge_t *existing_sgm,
				    const edge_t *existing_new_sgm);
static void purge_model(vcn_model_t *model);
static void set_intersection_holes(vcn_model_t *model,
				   const vcn_model_t *const model1,
				   const vcn_model_t *const model2);
static double* get_centroids_of_model_subareas(const vcn_model_t *model,
						uint32_t *N_centroids);
static uint32_t mask_intersection_holes(const vcn_model_t *model1,
					const vcn_model_t *model2,
					uint32_t N_centroids,
					const double *centroids,
					char *mask_centroids);
static void set_new_holes_to_model(uint32_t N_new_holes,
				   uint32_t N_centroids,
				   const double *centroids,
				   const char *mask_centroids,
				   vcn_model_t *model);
static void delete_isolated_elements(vcn_model_t *model);
static void delete_isolated_internal_vtx(vcn_model_t *model);
static void set_difference_holes(vcn_model_t *model,
				 const vcn_model_t *const model1,
				 const vcn_model_t *const model2);
static uint32_t mask_difference_holes(const vcn_model_t *model1,
				      const vcn_model_t *model2,
				      uint32_t N_centroids,
				      const double *centroids,
				      char *mask_centroids);
static void set_substraction_holes(vcn_model_t *model,
				   const vcn_model_t *const model2);
static uint32_t mask_substraction_holes(const vcn_model_t *model2,
					uint32_t N_centroids,
					const double *centroids,
					char *mask_centroids);

void vcn_model_get_combination(vcn_model_t *model,
			       const vcn_model_t *const input_model1,
			       const vcn_model_t *const input_model2)
{
	vcn_model_clear(model);

	vcn_model_t* model1 = (vcn_model_t*) input_model1;
	vcn_model_t* model2 = (vcn_model_t*) input_model2;

	scaling_data scaling;
	get_scaling_data(model1, model2, &scaling);
	scale_and_displace(model1, &scaling);
	scale_and_displace(model2, &scaling);

	/* Insert segments and vertices */
	nb_container_t* ht_vtx = alloca(nb_container_get_memsize(NB_HASH));
	nb_container_init(ht_vtx, NB_HASH);
	nb_container_set_key_generator(ht_vtx, vtx_hash_key);
	nb_container_set_comparer(ht_vtx, vtx_compare);
	nb_container_set_destroyer(ht_vtx, vtx_destroy);

	nb_container_t* avl_sgm = alloca(nb_container_get_memsize(NB_SORTED));
	nb_container_init(avl_sgm, NB_SORTED);
	nb_container_set_comparer(avl_sgm, edge_compare);
	nb_container_set_destroyer(avl_sgm, edge_destroy);

	get_combined_elements(model1, model2, ht_vtx, avl_sgm);

	if (nb_container_get_length(avl_sgm) > 2) {
		set_as_initial_vtx_in_edges(avl_sgm);
		delete_unused_vertices(ht_vtx);

		set_combined_model(model, ht_vtx, avl_sgm,
				   model1, model2);
		rescale_model(model, &scaling);
	}

	rescale_model(model1, &scaling);
	rescale_model(model2, &scaling);
	nb_container_finish(ht_vtx);
	nb_container_finish(avl_sgm);
}

static void get_combined_elements(const vcn_model_t *model1,
				  const vcn_model_t *model2,
				  nb_container_t *ht_vtx,
				  nb_container_t *avl_sgm)
{	   
	nb_container_t *sgm_intersect = 
		alloca(nb_container_get_memsize(NB_SORTED));
	nb_container_init(sgm_intersect, NB_SORTED);
	nb_container_set_comparer(sgm_intersect, ipack_comparer);

	search_intersections(sgm_intersect, model1, model2,
			     ht_vtx, avl_sgm);
	do {
		process_segment_intersections(avl_sgm, sgm_intersect, ht_vtx);
		remove_short_segments(avl_sgm, ht_vtx, MIN_LENGTH_X_SGM);
		search_intersections_in_edges(sgm_intersect, avl_sgm);
	} while (nb_container_is_not_empty(sgm_intersect));

	nb_container_finish(sgm_intersect);
}

static void search_intersections_in_edges(const nb_container_t *edges,
					  nb_container_t *intersections)
{  
	nb_iterator_t* iter = alloca(nb_iterator_get_memsize());
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, edges);
	while (nb_iterator_has_more(iter)) {
		const edge_t* edge1 = nb_iterator_get_next(iter);
		nb_iterator_t* subiter = alloca(nb_iterator_get_memsize());
		nb_iterator_copy(subiter, iter);
		while (nb_iterator_has_more(subiter)) {
			const edge_t* edge2 = nb_iterator_get_next(subiter);
			ipack_t *ipack = get_intersection_pack(edge1, edge2);
			if (NULL != ipack)
				nb_container_insert(intersections, ipack);
		}
		nb_iterator_finish(subiter);
	}
	nb_iterator_finish(iter);
}

static void set_combined_model(vcn_model_t *model,
			       nb_container_t *ht_vtx,
			       nb_container_t *avl_sgm,
			       const vcn_model_t *model1,
			       const vcn_model_t *model2)
{
	model->N = nb_container_get_length(ht_vtx);
	nb_model_alloc_vertices(model);

	model->M = nb_container_get_length(avl_sgm);
	nb_model_alloc_edges(model);

	set_vertices(model, ht_vtx);
	set_segments(model, avl_sgm);
	set_holes(model1, model2, model);
}

static void set_vertices(vcn_model_t *model,
			 nb_container_t *ht_vtx)
{
	nb_iterator_t* iter = alloca(nb_iterator_get_memsize());
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, ht_vtx);

	uint32_t vtx_id = 0;
	while (nb_iterator_has_more(iter)) {
		vtx_t *vtx = (vtx_t*) nb_iterator_get_next(iter);
		vtx_set_id(vtx, vtx_id);
		memcpy(GET_PVTX(model, vtx_id), vtx->x,
		       2 * sizeof(*(vtx->x)));
		vtx_id++;
	}
	nb_iterator_finish(iter);
}

static void set_segments(vcn_model_t *model,
			 nb_container_t *avl_sgm)
{
	nb_iterator_t* iter = alloca(nb_iterator_get_memsize());
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, avl_sgm);

	uint32_t sgm_counter = 0;
	while (nb_iterator_has_more(iter)) {
		const edge_t* edge = nb_iterator_get_next(iter);
		model->edge[sgm_counter * 2] = vtx_get_id(edge->v1);
		model->edge[sgm_counter*2+1] = vtx_get_id(edge->v2);
		sgm_counter++;		
	}
	nb_iterator_finish(iter);
}

static void set_holes(const vcn_model_t *model1,
		      const vcn_model_t *model2,
		      vcn_model_t *model)
{
	nb_mesh_t* mesh = alloca(nb_mesh_get_memsize());
	nb_mesh_init(mesh);

	nb_mesh_get_simplest_from_model(mesh, model);

	uint32_t N_centroids;
	double *centroids = 
		nb_mesh_get_centroids_of_subareas(mesh, &N_centroids);
	nb_mesh_finish(mesh);

	if (1 < N_centroids) {
		char* mask_centroids = nb_soft_allocate_mem(N_centroids);
		model->H = mask_true_centroids(model1, model2,
					       mask_centroids,
					       centroids,
					       N_centroids);

		if (model->H > 0) {
			nb_model_alloc_holes(model);
			model->H = 0;
			for (uint32_t i = 0; i < N_centroids; i++) {
				if (mask_centroids[i]) {
					memcpy(&(model->holes[model->H * 2]),
					       &(centroids[i * 2]),
					       2 * sizeof(*(model->holes)));
					model->H += 1;
				}
			}
		}
		nb_soft_free_mem(N_centroids, mask_centroids);
	}
	if (0 < N_centroids)
		free(centroids);
}

static uint32_t mask_true_centroids(const vcn_model_t *model1,
				    const vcn_model_t *model2,
				    char *mask, const double *centroids,
				    uint32_t N_centroids)
{
	memset(mask, 0, N_centroids);

	nb_mesh_t* mesh1 = alloca(nb_mesh_get_memsize());
	nb_mesh_init(mesh1);
	nb_mesh_get_simplest_from_model(mesh1, model1);

	nb_mesh_t* mesh2 = alloca(nb_mesh_get_memsize());
	nb_mesh_init(mesh2);
	nb_mesh_get_simplest_from_model(mesh2, model2);

	uint32_t N = 0;
	for (uint32_t i = 0; i < N_centroids; i++) {
		if (!nb_mesh_is_vtx_inside(mesh1, &(centroids[i*2]))) {
			if (!nb_mesh_is_vtx_inside(mesh2, &(centroids[i*2]))) {
				mask[i] = 1;
				N += 1;
			}
		}
	}
	nb_mesh_finish(mesh1);
	nb_mesh_finish(mesh2);
	return N;
}

static void rescale_model(vcn_model_t *model, const scaling_data *scaling)
{
	for (uint32_t i = 0; i < model->N; i++) {
		model->vertex[i * 2] = 
			model->vertex[i * 2] / scaling->scale + scaling->xdisp;
		model->vertex[i*2+1] =
			model->vertex[i*2+1] / scaling->scale + scaling->ydisp;
	}
  
	for (uint32_t i = 0; i < model->H; i++) {
		model->holes[i * 2] =
			model->holes[i * 2] / scaling->scale + scaling->xdisp;
		model->holes[i*2+1] =
			model->holes[i*2+1] / scaling->scale + scaling->ydisp;
	}
}

static void* ipack_create(void)
{
	return nb_allocate_zero_mem(sizeof(ipack_t));
}

static void ipack_destroy(void *ipack_ptr)
{
	free(ipack_ptr);
}

static int8_t ipack_comparer(const void *ipack1_ptr,
			     const void *ipack2_ptr)
{
	const ipack_t* ipack1 = ipack1_ptr;
	const ipack_t* ipack2 = ipack2_ptr;

	int8_t out;
	if (ipack1->status == ipack2->status) {
		if (ipack1 < ipack2)
			out = -1;
		else if (ipack1 > ipack2)
			out = 1;
		else
			out = 0;
	} else {
		if (NB_INTERSECTED == ipack1->status ||
		    NB_NOT_INTERSECTED == ipack2->status) {
			out = -1;
		} else if (NB_INTERSECTED == ipack2->status ||
			   NB_NOT_INTERSECTED == ipack1->status) {
			out = 1;
		} else if (NB_PARALLEL == ipack2->status) {
			out = -1;
		} else if (NB_PARALLEL == ipack1->status) {
			out = 1;
		} else {
			if (ipack1->status < ipack2->status)
				out = -1;
			else
				out = 1;
		}
	}
	return out;
}

static void get_scaling_data(const vcn_model_t *const model1,
			     const vcn_model_t *const model2,
			     scaling_data *scaling)
{
	double xmin = model1->vertex[0];
	double xmax = model1->vertex[0];
	double ymin = model1->vertex[1];
	double ymax = model1->vertex[1];
	for (uint32_t i = 1; i < model1->N; i++) {
		if (model1->vertex[i * 2] > xmax)
			xmax = model1->vertex[i * 2];
		else if (model1->vertex[i * 2] < xmin)
			xmin = model1->vertex[i * 2];

		if (model1->vertex[i*2+1] > ymax)
			ymax = model1->vertex[i*2+1];
		else if (model1->vertex[i*2+1] < ymin)
			ymin = model1->vertex[i*2+1];
	}
	for (uint32_t i = 0; i < model2->N; i++) {
		if (model2->vertex[i * 2] > xmax)
			xmax = model2->vertex[i * 2];
		else if (model2->vertex[i * 2] < xmin)
			xmin = model2->vertex[i * 2];
		if (model2->vertex[i*2+1] > ymax)
			ymax = model2->vertex[i*2+1];
		else if (model2->vertex[i*2+1] < ymin)
			ymin = model2->vertex[i*2+1];
	}

	scaling->xdisp = (xmin + xmax) / 2.0;
	scaling->ydisp = (ymin + ymax) / 2.0;
	scaling->scale = 100.0 / MAX(xmax-xmin, ymax-ymin);
}

static void scale_and_displace(vcn_model_t *model, 
			       const scaling_data *scaling)
{
	for (uint32_t i = 0; i < model->N; i++) {
		model->vertex[i * 2] = scaling->scale *
			(model->vertex[i * 2] - scaling->xdisp);
		model->vertex[i*2+1] = scaling->scale *
			(model->vertex[i*2+1] - scaling->ydisp);
	}
  
	for (uint32_t i = 0; i < model->H; i++) {
		model->holes[i * 2] = scaling->scale * 
			(model->holes[i * 2] - scaling->xdisp);
		model->holes[i*2+1] = scaling->scale *
			(model->holes[i*2+1] - scaling->ydisp);
	}
}

static void search_intersections(nb_container_t *intersections,
				 const vcn_model_t *model1,
				 const vcn_model_t *model2,
				 nb_container_t *vertices,
				 nb_container_t *edges)
{
	uint32_t memsize = (model1->M + model2->M) * sizeof(edge_t*);
	char *memblock = nb_soft_allocate_mem(memsize);
	memset(memblock, 0, memsize);

	edge_t **edges1 = (void*) memblock;
	edge_t **edges2 = (void*) (memblock + model1->M * sizeof(edge_t*));
	insert_edges_and_vtx(model1, vertices, edges, edges1, false);
	insert_edges_and_vtx(model2, vertices, edges, edges2, true);

	for (uint32_t i = 0; i < model1->M; i++) {
		for (uint32_t j = 0; j < model2->M; j++) {
			if (NULL != edges2[j]) {
				ipack_t* ipack = 
					get_intersection_pack(edges1[i],
							      edges2[j]);
				if (NULL != ipack)
					nb_container_insert(intersections,
							    ipack);
			}
		}
	}
	nb_soft_free_mem(memsize, memblock);
}

static void insert_edges_and_vtx(const vcn_model_t *model,
				 nb_container_t *vertices,
				 nb_container_t *edges,
				 edge_t **edges_array,
				 bool check_if_edges_exist)
{
	for (uint32_t i = 0; i < model->M; i++) {
		edge_t aux_edge;
		aux_edge.v1 = insert_vtx_if_not_exist(model, vertices,
						      GET_1_EDGE_VTX(model, i));
		aux_edge.v2 = insert_vtx_if_not_exist(model, vertices,
						      GET_2_EDGE_VTX(model, i));
		edge_set_length(&aux_edge);
		bool allow_to_insert_edge = true;
		if (check_if_edges_exist) {
			edge_t* existing = nb_container_exist(edges, &aux_edge);
			if (NULL != existing)
				allow_to_insert_edge = false;

		}
		if (allow_to_insert_edge) {
			edge_t *edge = edge_clone(&aux_edge);
			edges_array[i] = edge;
			nb_container_insert(edges, edge);
		} else {
			edges_array[i] = NULL;
		}
	}
}

static vtx_t* insert_vtx_if_not_exist(const vcn_model_t *const model,
				      nb_container_t *vertices,
				      uint32_t vtx_id)
{
	vtx_t vtx_to_find;
	memcpy(vtx_to_find.x, GET_PVTX(model, vtx_id), 
	       2 * sizeof(*(vtx_to_find.x)));

	vtx_t *vtx = nb_container_exist(vertices, &vtx_to_find);
	if (NULL == vtx) {
		vtx = vtx_create();
		memcpy(vtx->x, GET_PVTX(model, vtx_id), 2 * sizeof(*(vtx->x)));
		nb_container_insert(vertices, vtx);
	}
	return vtx;
}

static inline edge_t* exist_edge_and_insert_if_not(nb_container_t* edges, 
						   edge_t* edge)
{
	edge_t* existing_edge = nb_container_exist(edges, edge);
	if (NULL == existing_edge)
		nb_container_insert(edges, edge);
	return existing_edge;
}

static ipack_t* get_intersection_pack(const edge_t *sgm1,
				      const edge_t *sgm2)
{
	double intersection[2];
	nb_intersect_t status = vcn_utils2D_get_sgm_intersection(sgm1->v1->x,
								 sgm1->v2->x,
								 sgm2->v1->x,
								 sgm2->v2->x,
								 intersection);
	ipack_t *ipack = NULL;
	if (is_true_intersection(sgm1, sgm2, status)) {
		ipack = ipack_create();
		ipack->status = status;
		ipack->sgm1 = (edge_t*) sgm1;
		ipack->sgm2 = (edge_t*) sgm2;
		memcpy(ipack->intersection, intersection,
		       2 * sizeof(*intersection));
	}
	return ipack;
}

static bool is_true_intersection(const edge_t *sgm1, const edge_t *sgm2,
				 nb_intersect_t status)
{
	bool out = true;
	switch (status) {
	case NB_PARALLEL:
		if (if_collinear_sgm_are_not_intersected(sgm1, sgm2))
			out = false;
		break;
	case NB_INTERSECT_ON_A1:
		if (sgm1->v1 == sgm2->v1 || sgm1->v1 == sgm2->v2)
			out = false;
		break;
	case NB_INTERSECT_ON_A2:
		if (sgm1->v2 == sgm2->v1 || sgm1->v2 == sgm2->v2)
			out = false;
		break;
	case NB_INTERSECT_ON_B1:
		if (sgm2->v1 == sgm1->v1 || sgm2->v1 == sgm1->v2)
			out = false;
		break;
	case NB_INTERSECT_ON_B2:
		if (sgm2->v2 == sgm1->v1 || sgm2->v2 == sgm1->v2)
			out = false;
		break;
	}
	return out;
}

static bool if_collinear_sgm_are_not_intersected(const edge_t *sgm1,
						 const edge_t *sgm2)
{
	bool out = true;
	/* Check if are collinear */
	double orient = vcn_utils2D_orient(sgm1->v1->x,
					   sgm1->v2->x,
					   sgm2->v1->x);
	if (fabs(orient) < NB_GEOMETRIC_TOL) {
		collinear_dists cdists;
		calculate_collinear_dists(sgm1, sgm2, &cdists);

		if (collinear_are_intersected(&cdists))
			out = false;
	}
	return out;
}

static void process_segment_intersections(nb_container_t* avl_sgm,
					  nb_container_t* sgm_intersect,
					  nb_container_t* ht_vtx)
{
	while (nb_container_is_not_empty(sgm_intersect)) {
		ipack_t *ipack = nb_container_delete_first(sgm_intersect);
		process_ipack(ipack, avl_sgm, sgm_intersect, ht_vtx);
		ipack_destroy(ipack);
	}
}

static void process_ipack(ipack_t* ipack,
			  nb_container_t* avl_sgm,
			  nb_container_t* sgm_intersect,
			  nb_container_t* ht_vtx)
{
	switch (ipack->status) {
	case NB_INTERSECTED:
		process_intersected_sgm(ipack, avl_sgm,
					sgm_intersect, ht_vtx);
		break;
	case NB_PARALLEL:
		process_collinear_sgm(ipack, avl_sgm, sgm_intersect);
		break;
	case NB_INTERSECT_ON_A1:
		split_segment_by_vertex(avl_sgm, sgm_intersect, 
					ipack->sgm2, ipack->sgm1->v1);
		break;
    	case NB_INTERSECT_ON_A2:
		split_segment_by_vertex(avl_sgm, sgm_intersect,
					ipack->sgm2, ipack->sgm1->v2);
		break;
	case NB_INTERSECT_ON_B1:
		split_segment_by_vertex(avl_sgm, sgm_intersect,
					ipack->sgm1, ipack->sgm2->v1);
		break;
	case NB_INTERSECT_ON_B2:
		split_segment_by_vertex(avl_sgm, sgm_intersect,
					ipack->sgm1, ipack->sgm2->v2);
		break;
	}
}

static void process_intersected_sgm(ipack_t *ipack,
				    nb_container_t* avl_sgm,
				    nb_container_t* sgm_intersect,
				    nb_container_t* ht_vtx)
{
	/* Segments intersecting */
	vtx_t* vtx = alloca(vtx_get_memsize());
	vtx_init(vtx);
	memcpy(vtx->x, ipack->intersection, 2 * sizeof(*(vtx->x)));
      	vtx_t* aux_vtx = nb_container_exist(ht_vtx, vtx);
	vtx_finish(vtx);
	if (NULL == aux_vtx) {
		vtx = vtx_create();
		memcpy(vtx->x, ipack->intersection, 2 * sizeof(*(vtx->x)));
		nb_container_insert(ht_vtx, vtx);
	} else {
		vtx = aux_vtx;
	}

	/* Remove segments */
	nb_container_delete(avl_sgm, ipack->sgm1);
	nb_container_delete(avl_sgm, ipack->sgm2);

	split_segment_by_vertex(avl_sgm, sgm_intersect, ipack->sgm1, vtx);
	split_segment_by_vertex(avl_sgm, sgm_intersect, ipack->sgm2, vtx);
}

static void process_collinear_sgm(ipack_t *ipack,
				 nb_container_t* avl_sgm,
				 nb_container_t* sgm_intersect)
{
	collinear_dists cdists;
	calculate_collinear_dists(ipack->sgm1, ipack->sgm2, &cdists);
	if (collinear_are_contained(&cdists))
		process_collinear_contained(ipack, &cdists, avl_sgm,
					    sgm_intersect);
	else
		process_collinear_intersected(ipack, &cdists,
					      avl_sgm, sgm_intersect);
}

static void calculate_collinear_dists(const edge_t *sgm1,
				      const edge_t *sgm2,
				      collinear_dists *cdists)
{
	cdists->length1 = edge_get_length(sgm1);
	cdists->length2 = edge_get_length(sgm2);
	cdists->dist_1A_2A = vcn_utils2D_get_dist(sgm1->v1->x,
						  sgm2->v1->x);
	cdists->dist_1A_2B = vcn_utils2D_get_dist(sgm1->v1->x,
						  sgm2->v2->x);
	cdists->dist_1B_2A = vcn_utils2D_get_dist(sgm1->v2->x,
						  sgm2->v1->x);
	cdists->dist_1B_2B = vcn_utils2D_get_dist(sgm1->v2->x,
						  sgm2->v2->x);
	cdists->max_dist = MAX(cdists->dist_1B_2A, cdists->dist_1B_2B);
	cdists->max_dist = MAX(cdists->dist_1A_2B, cdists->max_dist);
	cdists->max_dist = MAX(cdists->dist_1A_2A, cdists->max_dist);
}

static inline bool collinear_are_contained(collinear_dists *cdists)
{
	double max_length = MAX(cdists->length1, cdists->length2);
	return (cdists->max_dist - max_length < -NB_GEOMETRIC_TOL);

}

static inline bool collinear_are_intersected(collinear_dists *cdists)
{
	double lengths_sum = cdists->length1 + cdists->length2;
	return (cdists->max_dist - lengths_sum < -NB_GEOMETRIC_TOL);

}

static void process_collinear_intersected(ipack_t *ipack,
					  collinear_dists *cdists,
					  nb_container_t* avl_sgm,
					  nb_container_t* sgm_intersect)
{
	/* Intersected by one side */
	nb_container_t* sgm_intersect_aux = 
		alloca(nb_container_get_memsize(NB_QUEUE));
	nb_container_init(sgm_intersect_aux, NB_QUEUE);

	nb_container_delete(avl_sgm, ipack->sgm1);
	nb_container_delete(avl_sgm, ipack->sgm2);

	collinear_merge_into_sgm1(ipack->sgm1, ipack->sgm2, cdists);
	edge_set_length(ipack->sgm1);

	edge_t* existing_sgm1 = 
		exist_edge_and_insert_if_not(avl_sgm, ipack->sgm1);

	update_ipacks_sgm1(ipack, sgm_intersect, sgm_intersect_aux,
			   existing_sgm1);	 
	update_ipacks_sgm2(ipack, sgm_intersect, sgm_intersect_aux,
			   existing_sgm1);
	    
	move_out_ipacks_from_aux(sgm_intersect, sgm_intersect_aux);

	if (NULL != existing_sgm1)
		edge_destroy(ipack->sgm1);
	edge_destroy(ipack->sgm2);
	nb_container_finish(sgm_intersect_aux);
}

static void collinear_merge_into_sgm1(edge_t *sgm1, const edge_t *sgm2, 
				      const collinear_dists *cdists)
{
	if (fabs(cdists->dist_1A_2A - cdists->max_dist) < 1e-8)
		sgm1->v2 = sgm2->v1;
	else if (fabs(cdists->dist_1A_2B - cdists->max_dist) < 1e-8)
		sgm1->v2 = sgm2->v2;
	else if (fabs(cdists->dist_1B_2A - cdists->max_dist) < 1e-8)
		sgm1->v1 = sgm2->v1;
	else
		sgm1->v1 = sgm2->v2;
}

static void update_ipacks_sgm1(ipack_t *ipack,
			       nb_container_t *sgm_intersect,
			       nb_container_t *sgm_intersect_aux,
			       const edge_t *existing_sgm1)
{
	ipack_t* jpack = search_ipack_with_sgm(sgm_intersect, ipack->sgm1);
	while (NULL != jpack) {
		nb_container_delete(sgm_intersect, jpack);
		edge_t* sgm_aux = ipack_get_opp_sgm(jpack, ipack->sgm1);
		ipack_destroy(jpack);

		if (NULL == existing_sgm1) {
			jpack = get_intersection_pack(ipack->sgm1, sgm_aux);
			if (NULL != jpack)
				nb_container_insert(sgm_intersect_aux, jpack);
		}

	      	jpack = search_ipack_with_sgm(sgm_intersect, ipack->sgm1);
	}
}

static edge_t *ipack_get_opp_sgm(const ipack_t *ipack,
					 const edge_t *sgm)
{
	edge_t* sgm_aux = ipack->sgm1;
	if (sgm == sgm_aux)
		sgm_aux = ipack->sgm2;
	return sgm_aux;

}

static void update_ipacks_sgm2(ipack_t *ipack,
			       nb_container_t *sgm_intersect,
			       nb_container_t *sgm_intersect_aux,
			       const edge_t *existing_sgm1)
{
	ipack_t *jpack = search_ipack_with_sgm(sgm_intersect, ipack->sgm2);
	while (NULL != jpack) {
		nb_container_delete(sgm_intersect, jpack);
		edge_t* sgm_aux = ipack_get_opp_sgm(jpack, ipack->sgm2);
		ipack_destroy(jpack);
	      
		if (NULL == existing_sgm1) {
			jpack = search_ipack_with_both_sgm(sgm_intersect_aux,
							   ipack->sgm1,
							   sgm_aux);
			if (NULL == jpack) {
				jpack = get_intersection_pack(ipack->sgm1,
							      sgm_aux);

				if (NULL != jpack)
					nb_container_insert(sgm_intersect_aux,
							    jpack);
			}
		}
	      
		jpack = search_ipack_with_sgm(sgm_intersect, ipack->sgm2);
	}

}

static void move_out_ipacks_from_aux(nb_container_t *sgm_intersect,
				     nb_container_t *sgm_intersect_aux)
{
	while (nb_container_is_not_empty(sgm_intersect_aux)) {
		ipack_t* pack = nb_container_delete_first(sgm_intersect_aux);
		nb_container_insert(sgm_intersect, pack);
	}
}

static void process_collinear_contained(ipack_t *ipack,
					collinear_dists *cdists,
					nb_container_t* avl_sgm,
					nb_container_t* sgm_intersect)
{
	edge_t* small_sgm = get_smallest_sgm(ipack, cdists);

	nb_container_delete(avl_sgm, small_sgm);

	ipack_t* jpack = search_ipack_with_sgm(sgm_intersect, small_sgm);
	while (NULL != jpack) {
		nb_container_delete(sgm_intersect, jpack);
		ipack_destroy(jpack);
		jpack = search_ipack_with_sgm(sgm_intersect, small_sgm);
	}
	edge_destroy(small_sgm);
}

static edge_t *get_smallest_sgm(ipack_t *ipack, collinear_dists *cdists)
{
	edge_t* small_sgm;
	if (cdists->length1 > cdists->length2)
		small_sgm = ipack->sgm2;  
	else
		small_sgm = ipack->sgm1;
	return small_sgm;
}

static void remove_short_segments(nb_container_t* avl_sgm,
				  nb_container_t* ht_vtx,
				  double min_length_x_sgm)
{
	nb_container_t *short_sgm = alloca(nb_container_get_memsize(NB_QUEUE));
	nb_container_init(short_sgm, NB_QUEUE);
	get_short_sgm(avl_sgm, short_sgm, min_length_x_sgm);

	while (nb_container_is_not_empty(short_sgm)) {
		edge_t* sgm = nb_container_delete_first(short_sgm);
		nb_container_delete(avl_sgm, sgm);

		collapse_sgm_into_v1(avl_sgm, ht_vtx, sgm);
		/* Free memory */
		vtx_destroy(sgm->v2);
		edge_destroy(sgm);
	}
	nb_container_finish(short_sgm);
}

static void get_short_sgm(const nb_container_t *avl_sgm,
			  nb_container_t *short_sgm,
			  double min_length_x_sgm)
{
	nb_iterator_t *iter = alloca(nb_iterator_get_memsize());
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, avl_sgm);
	while (nb_iterator_has_more(iter)) {
		edge_t* sgm = (edge_t*) nb_iterator_get_next(iter);
		double length = edge_get_length(sgm);
		if (length < MAX(min_length_x_sgm, NB_GEOMETRIC_TOL))
			nb_container_insert(short_sgm, sgm);
	}
	nb_iterator_finish(iter);
}

static void collapse_sgm_into_v1(nb_container_t *avl_sgm,
				 nb_container_t *ht_vtx, edge_t *sgm)
{
	nb_container_delete(ht_vtx, sgm->v1);
	nb_container_delete(ht_vtx, sgm->v2);

	nb_container_t* list_sgm = alloca(nb_container_get_memsize(NB_QUEUE));
	nb_container_init(list_sgm, NB_QUEUE);

	get_connected_segments(avl_sgm, sgm, list_sgm);

	link_connected_segments_to_v1(sgm, avl_sgm, list_sgm);
	v1_to_centroid_if_connections(list_sgm, ht_vtx, sgm);

	update_segments_connected(avl_sgm, list_sgm);

	nb_container_finish(list_sgm);
}

static void link_connected_segments_to_v1(const edge_t *sgm,
					  nb_container_t *avl_sgm,
					  nb_container_t *list_sgm)
{
	nb_iterator_t* iter = alloca(nb_iterator_get_memsize());
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, list_sgm);
	while (nb_iterator_has_more(iter)) {
		edge_t* subsgm = (edge_t*) nb_iterator_get_next(iter);

		nb_container_delete(avl_sgm, subsgm);

		if (subsgm->v1 == sgm->v2)
			subsgm->v1 = sgm->v1;
		else if (subsgm->v2 == sgm->v2)
			subsgm->v2 = sgm->v1;
	}
	nb_iterator_finish(iter);

}

static void get_connected_segments(const nb_container_t *avl_sgm,
				   const edge_t *sgm,
				   nb_container_t *list_sgm)
{
	nb_iterator_t* iter = alloca(nb_iterator_get_memsize());
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, avl_sgm);
	while (nb_iterator_has_more(iter)) {
		const edge_t* subsgm = nb_iterator_get_next(iter);
		if (subsgm->v1 == sgm->v2 || subsgm->v2 == sgm->v2 ||
		    subsgm->v1 == sgm->v1 || subsgm->v2 == sgm->v1)
			nb_container_insert(list_sgm, subsgm);
	}
	nb_iterator_finish(iter);
}

static void v1_to_centroid_if_connections(const nb_container_t *list_sgm,
					  nb_container_t *ht_vtx,
					  edge_t *sgm)
{
	if (nb_container_is_not_empty(list_sgm)) {
		sgm->v1->x[0] = (sgm->v1->x[0] + sgm->v2->x[0]) * 0.5;
		sgm->v1->x[1] = (sgm->v1->x[1] + sgm->v2->x[1]) * 0.5;
		nb_container_insert(ht_vtx, sgm->v1);
	} else {
		vtx_destroy(sgm->v1);
	}
}

static void update_segments_connected(nb_container_t *avl_sgm,
				      nb_container_t *list_sgm)
{
	while (nb_container_is_not_empty(list_sgm)) {
		edge_t* subsgm = nb_container_delete_first(list_sgm);
		/* Update length */
		edge_set_length(subsgm);

		/* Insert if it does not exist */
		edge_t* existing_sgm =
			exist_edge_and_insert_if_not(avl_sgm, subsgm);
		if (NULL != existing_sgm)
			edge_destroy(subsgm);
	}
}

static void set_as_initial_vtx_in_edges(nb_container_t *edges)
{
	nb_iterator_t *iter = alloca(nb_iterator_get_memsize());
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, edges);
	while (nb_iterator_has_more(iter)) {
		const edge_t *edge = nb_iterator_get_next(iter);
		vtx_set_as_initial(edge->v1);
		vtx_set_as_initial(edge->v2);
	}
	nb_iterator_finish(iter);
}

static void delete_unused_vertices(nb_container_t *vertices)
{
	nb_container_t *to_delete = alloca(nb_container_get_memsize(NB_QUEUE));
	nb_container_init(to_delete, NB_QUEUE);
	nb_iterator_t *iter = alloca(nb_iterator_get_memsize());
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, vertices);
	while (nb_iterator_has_more(iter)) {
		const vtx_t* vtx = nb_iterator_get_next(iter);
		if (vtx_is_not_initial(vtx))
			nb_container_insert(to_delete, vtx);
	}
	nb_iterator_finish(iter);

	while (nb_container_is_not_empty(to_delete)) {
		vtx_t *vtx = nb_container_delete_first(to_delete);
		nb_container_delete(vertices, vtx);
		vtx_destroy(vtx);		
	}
	nb_container_finish(to_delete);
}

static ipack_t* search_ipack_with_sgm
				(const nb_container_t *const packs,
				 const edge_t *const edge)
{
	nb_iterator_t* iter = alloca(nb_iterator_get_memsize());
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, packs);
	ipack_t *out_ipack = NULL;
	while (nb_iterator_has_more(iter)) {
		ipack_t *ipack = (ipack_t*) nb_iterator_get_next(iter);
		if (edge == ipack->sgm1 || edge == ipack->sgm2) {
			out_ipack = ipack;
			break;
		}
	}
	nb_iterator_finish(iter);
	return out_ipack;
}

static ipack_t* search_ipack_with_both_sgm(const nb_container_t *const packs,
					   const edge_t *const edge1,
					   const edge_t *const edge2)
{
	nb_iterator_t* iter = alloca(nb_iterator_get_memsize());
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, packs);
	ipack_t *out_ipack = NULL;
	while (nb_iterator_has_more(iter)) {
		ipack_t *ipack = (ipack_t*) nb_iterator_get_next(iter);
		if ((edge1 == ipack->sgm1 && edge2 == ipack->sgm2) ||
		    (edge2 == ipack->sgm1 && edge1 == ipack->sgm2)) {
			out_ipack = ipack;
			break;
		}
	}
	nb_iterator_finish(iter);
	return out_ipack;
}

static void split_segment_by_vertex(nb_container_t* avl_sgm,
				    nb_container_t* sgm_intersect,
				    edge_t* sgm, vtx_t* vtx)
{
	nb_container_delete(avl_sgm, sgm);
	edge_t* new_sgm = edge_create();
	
	new_sgm->v1 = vtx;
	new_sgm->v2 = sgm->v2;
	sgm->v2 = vtx;

	/* Update lengths */
	edge_set_length(sgm);
	edge_set_length(new_sgm);

	/* Insert into the AVL */
	edge_t* existing_sgm = exist_edge_and_insert_if_not(avl_sgm, sgm);

	edge_t* existing_new_sgm =
		exist_edge_and_insert_if_not(avl_sgm, new_sgm);
	if (NULL != existing_new_sgm)
		edge_destroy(new_sgm);
	
	nb_container_t* sgm_intersect_aux =
		alloca(nb_container_get_memsize(NB_QUEUE));
	nb_container_init(sgm_intersect_aux, NB_QUEUE);
  
	set_new_segments_ipacks(sgm, new_sgm, 
				sgm_intersect, sgm_intersect_aux,
				existing_sgm, existing_new_sgm);	 

	move_out_ipacks_from_aux(sgm_intersect,  sgm_intersect_aux);

	nb_container_finish(sgm_intersect_aux);

	if (NULL != existing_sgm)
		edge_destroy(sgm);
}

static void set_new_segments_ipacks(const edge_t *sgm,
				    const edge_t *new_sgm,
				    nb_container_t *sgm_intersect,
				    nb_container_t *sgm_intersect_aux,
				    const edge_t *existing_sgm,
				    const edge_t *existing_new_sgm)
{
	ipack_t* ipack = search_ipack_with_sgm(sgm_intersect, sgm);
	while (NULL != ipack) {
		edge_t* sgm_aux = ipack_get_opp_sgm(ipack, sgm);

		if (NULL == existing_sgm) {
			ipack_t *new_pack = get_intersection_pack(sgm, sgm_aux);
			if (NULL != new_pack)
				nb_container_insert(sgm_intersect_aux, new_pack);
		}

		if (NULL == existing_new_sgm) {
			ipack_t *new_pack = get_intersection_pack(new_sgm,
								  sgm_aux);
			if (NULL != new_pack)
				nb_container_insert(sgm_intersect_aux, new_pack);
		}
		nb_container_delete(sgm_intersect, ipack);

		ipack_destroy(ipack);

		ipack = search_ipack_with_sgm(sgm_intersect, sgm);
	}
}

void vcn_model_get_intersection(vcn_model_t *model,
				const vcn_model_t *const model1,
				const vcn_model_t *const model2)
{
	vcn_model_get_combination(model, model1, model2);
	set_intersection_holes(model, model1, model2);
	purge_model(model);

}

static void purge_model(vcn_model_t *model)
{
	delete_isolated_elements(model);
	if (0 < model->N)
		delete_isolated_internal_vtx(model);
}

static void set_intersection_holes(vcn_model_t *model,
				   const vcn_model_t *const model1,
				   const vcn_model_t *const model2)
{
	uint32_t N_centroids;
	double *centroids = 
		get_centroids_of_model_subareas(model, &N_centroids);

	if (0 < N_centroids) {
		char* mask_centroids = nb_soft_allocate_mem(N_centroids);

		uint32_t N_new_holes =
			mask_intersection_holes(model1, model2, N_centroids,
						centroids, mask_centroids);
	
		set_new_holes_to_model(N_new_holes, N_centroids,
				       centroids, mask_centroids, model);
		nb_soft_free_mem(N_centroids, mask_centroids);
		free(centroids);
	}
}

static double* get_centroids_of_model_subareas(const vcn_model_t *model,
						uint32_t *N_centroids)
{
	nb_mesh_t* mesh = alloca(nb_mesh_get_memsize());
	nb_mesh_init(mesh);

	nb_mesh_get_simplest_from_model(mesh, model);
	double* centroids = 
		nb_mesh_get_centroids_of_subareas(mesh, N_centroids);
	nb_mesh_finish(mesh);
	return centroids;
}

static uint32_t mask_intersection_holes(const vcn_model_t *model1,
					const vcn_model_t *model2,
					uint32_t N_centroids,
					const double *centroids,
					char *mask_centroids)
{
	uint32_t N_new_holes = 0;
	nb_mesh_t* mesh1 = alloca(nb_mesh_get_memsize());
	nb_mesh_init(mesh1);
	nb_mesh_get_simplest_from_model(mesh1, model1);

	nb_mesh_t* mesh2 = alloca(nb_mesh_get_memsize());
	nb_mesh_init(mesh2);
	nb_mesh_get_simplest_from_model(mesh2, model2);
 
	for (uint32_t i = 0; i < N_centroids; i++) {
		if (!nb_mesh_is_vtx_inside(mesh1, &(centroids[i * 2])) ||
		    !nb_mesh_is_vtx_inside(mesh2, &(centroids[i * 2]))) {
			mask_centroids[i] = 1;
			N_new_holes += 1;
		} else {
			mask_centroids[i] = 0;
		}
	}
	nb_mesh_finish(mesh1);
	nb_mesh_finish(mesh2);
	return N_new_holes;
}

static void set_new_holes_to_model(uint32_t N_new_holes,
				   uint32_t N_centroids,
				   const double *centroids,
				   const char *mask_centroids,
				   vcn_model_t *model)
{
	if (0 < N_new_holes) {
		uint32_t total_holes = model->H + N_new_holes;
		double* new_holes = nb_allocate_zero_mem(2 * total_holes *
							 sizeof(*new_holes));
		if (0 < model->H) {
			memcpy(new_holes, model->holes,
			       2 * model->H * sizeof(*new_holes));
			free(model->holes);
		}

		uint32_t i = model->H;
		for (uint32_t j = 0; j < N_centroids; j++) {
			if (1 == mask_centroids[j]) {
				memcpy(&(new_holes[i * 2]), 
				       &(centroids[j * 2]),
				       2 * sizeof(*new_holes));
				i += 1;
			}
		}

		model->H = total_holes;
		model->holes = new_holes;
	}
}

static void delete_isolated_elements(vcn_model_t *model)
{
	uint32_t mesh_memsize = nb_mesh_get_memsize();
	nb_mesh_t *mesh = alloca(mesh_memsize);
	nb_mesh_init(mesh);

	nb_mesh_get_simplest_from_model(mesh, model);

	nb_mesh_delete_isolated_segments(mesh);
	nb_mesh_delete_internal_input_segments(mesh);
	nb_mesh_delete_isolated_vertices(mesh);

	if (0 == nb_mesh_get_N_trg(mesh)) {
		vcn_model_clear(model);
	} else {
		uint32_t msh_memsize = nb_partition_get_memsize(NB_TRIAN);
		nb_partition_t* part = alloca(msh_memsize);
		nb_partition_init(part, NB_TRIAN);
		nb_partition_load_from_mesh(part, mesh);

		nb_partition_build_model(part, model);
		nb_partition_finish(part);
	}
	nb_mesh_finish(mesh);
}

static void delete_isolated_internal_vtx(vcn_model_t *model)
{
	uint32_t mask_size = model->N;
	char* mask = nb_soft_allocate_mem(mask_size);
	memset(mask, 0, mask_size);
	for(uint32_t i = 0; i < 2 * model->M; i++)
		mask[model->edge[i]] = 1;

	uint32_t N_vtx = 0;
	uint32_t perm_memsize = model->N * sizeof(uint32_t);
	uint32_t* perm = nb_soft_allocate_mem(perm_memsize);
	for (uint32_t i = 0; i < model->N; i++){
		if (1 == mask[i]) {
			perm[i] = N_vtx;
			N_vtx ++;
		}
	}
	if (0 == N_vtx) {
		model->N = 0;
		free(model->vertex);
		model->vertex = NULL;
	} else if (N_vtx < model->N) {
		double* vertices = nb_allocate_zero_mem(N_vtx * 2 *
							sizeof(*vertices));
		for (uint32_t i = 0; i < model->N; i++) {
			if (1 == mask[i]) {
				memcpy(&(vertices[perm[i] * 2]),
				       &(model->vertex[i * 2]), 
				       2 * sizeof(double));
			}
		}  
    
		for (uint32_t i = 0; i < 2 * model->M; i++) {
			model->edge[i] = perm[model->edge[i]];
		}
		model->N = N_vtx;
		free(model->vertex);
		model->vertex = vertices;
	}
	nb_soft_free_mem(mask_size, mask);
	nb_soft_free_mem(perm_memsize, perm);
}

void vcn_model_get_union(vcn_model_t *model,
			 const vcn_model_t *model1,
			 const vcn_model_t *model2)
{
	vcn_model_get_combination(model, model1, model2);
	purge_model(model);
}


void vcn_model_get_difference(vcn_model_t *model,
			      const vcn_model_t *const model1,
			      const vcn_model_t *const model2)
{
	vcn_model_get_combination(model, model1, model2);
	set_difference_holes(model, model1, model2);
	purge_model(model);

}

static void set_difference_holes(vcn_model_t *model,
				 const vcn_model_t *const model1,
				 const vcn_model_t *const model2)
{
	uint32_t N_centroids;
	double *centroids = 
		get_centroids_of_model_subareas(model, &N_centroids);

	if (0 < N_centroids) {
		char* mask_centroids = nb_soft_allocate_mem(N_centroids);

		uint32_t N_new_holes =
			mask_difference_holes(model1, model2, N_centroids,
					      centroids, mask_centroids);
	
		set_new_holes_to_model(N_new_holes, N_centroids,
				       centroids, mask_centroids, model);
		nb_soft_free_mem(N_centroids, mask_centroids);
		free(centroids);
	}
}

static uint32_t mask_difference_holes(const vcn_model_t *model1,
				      const vcn_model_t *model2,
				      uint32_t N_centroids,
				      const double *centroids,
				      char *mask_centroids)
{
	uint32_t N_new_holes = 0;
	nb_mesh_t* mesh1 = alloca(nb_mesh_get_memsize());
	nb_mesh_init(mesh1);
	nb_mesh_get_simplest_from_model(mesh1, model1);

	nb_mesh_t* mesh2 = alloca(nb_mesh_get_memsize());
	nb_mesh_init(mesh2);
	nb_mesh_get_simplest_from_model(mesh2, model2);

	for (uint32_t i = 0; i < N_centroids; i++) {
		if (nb_mesh_is_vtx_inside(mesh1, &(centroids[i * 2])) &&
		    nb_mesh_is_vtx_inside(mesh2, &(centroids[i * 2]))) {
			mask_centroids[i] = 1;
			N_new_holes += 1;
		} else {
			mask_centroids[i] = 0;
		}
	}
	nb_mesh_finish(mesh1);
	nb_mesh_finish(mesh2);
	return N_new_holes;
}

void vcn_model_get_substraction(vcn_model_t *model,
				const vcn_model_t *const model1,
				const vcn_model_t *const model2)
{
	vcn_model_get_combination(model, model1, model2);
	set_substraction_holes(model, model2);
	purge_model(model);
}

static void set_substraction_holes(vcn_model_t *model,
				   const vcn_model_t *const model2)
{
	uint32_t N_centroids;
	double *centroids = 
		get_centroids_of_model_subareas(model, &N_centroids);

	if (0 < N_centroids) {
		char* mask_centroids = nb_soft_allocate_mem(N_centroids);
	
		uint32_t N_new_holes =
			mask_substraction_holes(model2, N_centroids,
						centroids, mask_centroids);

		set_new_holes_to_model(N_new_holes, N_centroids,
				       centroids, mask_centroids, model);

		free(centroids);
		nb_soft_free_mem(N_centroids, mask_centroids);
	}
}

static uint32_t mask_substraction_holes(const vcn_model_t *model2,
					uint32_t N_centroids,
					const double *centroids,
					char *mask_centroids)
{
	uint32_t N_new_holes = 0;

	nb_mesh_t* mesh2 = alloca(nb_mesh_get_memsize());
	nb_mesh_init(mesh2);
	nb_mesh_get_simplest_from_model(mesh2, model2);
	for (uint32_t i = 0; i < N_centroids; i++) {
		if (nb_mesh_is_vtx_inside(mesh2, &(centroids[i * 2]))) {
			mask_centroids[i] = 1;
			N_new_holes += 1;
		} else {
			mask_centroids[i] = 0;
		}
	}
	nb_mesh_finish(mesh2);
	return N_new_holes;
}
