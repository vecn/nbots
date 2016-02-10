#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>
#include <math.h>

#include "nb/container_bot/container.h"
#include "nb/container_bot/iterator.h"
#include "nb/geometric_bot/utils2D.h"
#include "nb/geometric_bot/mesh/mesh2D.h"
#include "nb/geometric_bot/mesh/elements2D/triangles.h"
#include "nb/geometric_bot/mesh/modules2D/pruner.h"

#include "vtx.h"
#include "edge.h"
#include "../model2D_struct.h"
#include "nb/geometric_bot/model/model2D.h"
#include "nb/geometric_bot/model/modules2D/verifier.h"
#include "nb/geometric_bot/model/modules2D/blender.h"

#define GET_PVTX(model, i) (&((model)->vertex[(i)*2]))
#define GET_1_EDGE_VTX(model, i) ((model)->edge[(i) * 2])
#define GET_2_EDGE_VTX(model, i) ((model)->edge[(i)*2+1])
#define MAX(a, b) (((a)>(b))?(a):(b))

static double get_scale_and_disp(const vcn_model_t *const model1,
				 const vcn_model_t *const model2,
				 double *xdisp, double *ydisp);
static void scale_and_displace(vcn_model_t *model, double scale,
			       double xdisp, double ydisp);
static vcn_container_t* search_intersections(vcn_model_t *model1,
					     vcn_model_t *model2,
					     vcn_container_t *vertices,
					     vcn_container_t *edges);
static edge_t** insert_edges_and_vtx(vcn_model_t *model,
				      vcn_container_t *vertices,
				      vcn_container_t *edges,
				      bool check_if_edges_exist);
static vtx_t* insert_vtx_if_not_exist(const vcn_model_t *const model,
				      vcn_container_t *vertices,
				      uint32_t vtx_id);
static uint32_t intersection_pack_key(const void* const pack);
static edge_t* exist_edge_and_insert_if_not(vcn_container_t* edges,
					    edge_t* edge);
static void** get_intersection_pack(const edge_t *const sgm1,
				    const edge_t *const sgm2);
static void process_segment_intersections(vcn_container_t* avl_sgm,
					  vcn_container_t* sgm_intersect,
					  vcn_container_t* ht_vtx);
static void process_short_edges(vcn_container_t *vertices,
				vcn_container_t *edges,
				double min_length_x_edg);
static void remove_short_segments(vcn_container_t* avl_sgm,
				  vcn_container_t* ht_vtx,
				  double min_length_x_segment);
static vcn_container_t* search_intersections_in_edges(vcn_container_t *edges);
static void set_as_initial_vtx_in_edges(vcn_container_t *edges);
static void delete_unused_vertices(vcn_container_t *vertices);
static void** search_intersection_pack_with_sgm
				(const vcn_container_t *const packs,
				 const edge_t *const edge);
static void destroy_intersection_pack(void **pack);
static void** search_intersection_pack_with_both_segments
			(const vcn_container_t *const packs,
			 const edge_t *const edge1,
			 const edge_t *const edge2);
static void split_segment_by_vertex(vcn_container_t* avl_sgm,
				    vcn_container_t* sgm_intersect,
				    edge_t* sgm, vtx_t* vtx);

vcn_model_t* vcn_model_get_combination(const vcn_model_t *const input_model1,
				       const vcn_model_t *const input_model2,
				       double min_length_x_segment)
{
	vcn_model_t* model1 = vcn_model_clone(input_model1);
	vcn_model_t* model2 = vcn_model_clone(input_model2);

	/* Calculate displacement and scale to avoid floating imprecisions */
	double xdisp, ydisp;
	double scale = get_scale_and_disp(model1, model2, &xdisp, &ydisp);
	scale_and_displace(model1, scale, xdisp, ydisp);
	scale_and_displace(model2, scale, xdisp, ydisp);

	/* Insert segments and vertices */
	vcn_container_t* ht_vtx = vcn_container_create(NB_CONTAINER_HASH);
	vcn_container_set_key_generator(ht_vtx, vtx_hash_key);
	vcn_container_set_comparer(ht_vtx, vtx_are_equal);
	vcn_container_set_destroyer(ht_vtx, vtx_destroy);

	vcn_container_t* avl_sgm = vcn_container_create(NB_CONTAINER_SORTED);
	vcn_container_set_key_generator(avl_sgm, edge_key_by_length);
	vcn_container_set_comparer(avl_sgm, edge_are_equal);
	vcn_container_set_destroyer(avl_sgm, edge_destroy);
	   
	vcn_container_t* sgm_intersect = search_intersections(model1, model2,
							      ht_vtx, avl_sgm);

	process_segment_intersections(avl_sgm, sgm_intersect, ht_vtx);

	vcn_container_destroy(sgm_intersect);

	process_short_edges(ht_vtx, avl_sgm, min_length_x_segment * scale);

	/* Check that a volume is possible */
	if (vcn_container_get_length(avl_sgm) < 3) {
		vcn_model_destroy(model1);
		vcn_model_destroy(model2);
		vcn_container_destroy(avl_sgm);
		vcn_container_destroy(ht_vtx);
		return vcn_model_create();
	}
	
	set_as_initial_vtx_in_edges(avl_sgm);
	delete_unused_vertices(ht_vtx);

	/* Allocate new model */
	vcn_model_t* model = vcn_model_create();
	model->N = vcn_container_get_length(ht_vtx);
	model_alloc_vertices(model);
	model->M = vcn_container_get_length(avl_sgm);
	model_alloc_edges(model);

	/* Set vertices to the model */
	uint32_t vtx_id = 0;
	vcn_iterator_t *ht_iter = vcn_iterator_create();
	vcn_iterator_set_container(ht_iter, ht_vtx);
	while(vcn_iterator_has_more(ht_iter)){
		const vtx_t* vtx = vcn_iterator_get_next(ht_iter);
		vtx_set_id((vtx_t*)vtx, vtx_id);
		memcpy(GET_PVTX(model, vtx_id), vtx->x, 2 * sizeof(*(vtx->x)));
		vtx_id++;
	}
	vcn_iterator_destroy(ht_iter);

	/* Set segments to the model */
	uint32_t sgm_counter = 0;
	vcn_iterator_t *avl_iter = vcn_iterator_create();
	vcn_iterator_set_container(avl_iter, avl_sgm);
	while(vcn_iterator_has_more(avl_iter)){
		const edge_t* edge = vcn_iterator_get_next(avl_iter);
		model->edge[sgm_counter * 2] = vtx_get_id(edge->v1);
		model->edge[sgm_counter*2+1] = vtx_get_id(edge->v2);
		sgm_counter++;
	}
	vcn_iterator_destroy(avl_iter);

	/* Get holes */
	vcn_mesh_t* mesh = vcn_mesh_create();
	vcn_mesh_generate_from_model(mesh, model);

	uint32_t N_centroids;
	double* centroids = vcn_mesh_get_centroids_of_subareas(mesh, &N_centroids);
	vcn_mesh_destroy(mesh);

	if (N_centroids > 1) {
		char* mask_centroids = calloc(N_centroids, 1);
		
		vcn_mesh_t* mesh1 = vcn_mesh_create();
		vcn_mesh_generate_from_model(mesh1, model1);

		vcn_mesh_t* mesh2 = vcn_mesh_create();
		vcn_mesh_generate_from_model(mesh2, model2);

		model->H = 0;
		for (uint32_t i = 0; i < N_centroids; i++) {
			if (vcn_mesh_is_vtx_inside(mesh1, &(centroids[i * 2])))
				continue;
			if (vcn_mesh_is_vtx_inside(mesh2, &(centroids[i * 2])))
				continue;
			mask_centroids[i] = 1;
			model->H += 1;
		}
		vcn_mesh_destroy(mesh1);
		vcn_mesh_destroy(mesh2);

		if (model->H > 0) {
			model_alloc_holes(model);
			model->H = 0;
			for (uint32_t i=0; i < N_centroids; i++) {
				if (mask_centroids[i] == 0)
					continue;
				memcpy(&(model->holes[model->H * 2]),
				       &(centroids[i * 2]), 2 * sizeof(*(model->holes)));
				model->H += 1;
			}
		}

		/* Free memory */
		free(mask_centroids);
	}
	free(centroids);
	vcn_model_destroy(model1);
	vcn_model_destroy(model2);
  
	vcn_container_set_destroyer(ht_vtx, (void (*)(void*)) vtx_destroy);
	vcn_container_destroy(ht_vtx);
  
	while (vcn_container_is_not_empty(avl_sgm)) {
		edge_t* sgm = vcn_container_delete_first(avl_sgm);
		edge_destroy(sgm);
	}
	vcn_container_destroy(avl_sgm);

	/* Rescale coordinates (to avoid floating point error) */
	for (uint32_t i = 0; i < model->N; i++) {
		model->vertex[i * 2] = model->vertex[i * 2]/scale + xdisp;
		model->vertex[i*2+1] = model->vertex[i*2+1]/scale + ydisp;
	}
  
	for (uint32_t i = 0; i < model->H; i++) {
		model->holes[i * 2] = model->holes[i * 2]/scale + xdisp;
		model->holes[i*2+1] = model->holes[i*2+1]/scale + ydisp;
	}

	/* Return models combination */
	return model;
}

static double get_scale_and_disp(const vcn_model_t *const model1,
				 const vcn_model_t *const model2,
				 double *xdisp, double *ydisp)
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

	*xdisp = (xmin + xmax) / 2.0;
	*ydisp = (ymin + ymax) / 2.0;
	double scale = 100.0 / MAX(xmax-xmin, ymax-ymin);
	return scale;
}

static void scale_and_displace(vcn_model_t *model, double scale,
			       double xdisp, double ydisp)
{
	for(uint32_t i = 0; i < model->N; i++){
		model->vertex[i * 2] = scale * (model->vertex[i * 2] - xdisp);
		model->vertex[i*2+1] = scale * (model->vertex[i*2+1] - ydisp);
	}
  
	for(uint32_t i = 0; i < model->H; i++){
		model->holes[i * 2] = scale * (model->holes[i * 2] - xdisp);
		model->holes[i*2+1] = scale * (model->holes[i*2+1] - ydisp);
	}
}

static vcn_container_t* search_intersections(vcn_model_t *model1,
					     vcn_model_t *model2,
					     vcn_container_t *vertices,
					     vcn_container_t *edges)
{
	edge_t** edges1 = insert_edges_and_vtx(model1, vertices, edges, false);
	edge_t** edges2 = insert_edges_and_vtx(model2, vertices, edges, true);
	vcn_container_t *intersections =
		vcn_container_create(NB_CONTAINER_SORTED);
	vcn_container_set_key_generator(intersections, intersection_pack_key);
	for (uint32_t i = 0; i < model1->M; i++) {
		for (uint32_t j = 0; j < model2->M; j++) {
			if (NULL != edges2[j]) {
				void** intersection_pack = 
					get_intersection_pack(edges1[i],
							      edges2[j]);
				if (NULL != intersection_pack)
					vcn_container_insert(intersections,
							     intersection_pack);
			}
		}
	}
	free(edges1);
	free(edges2);
	return intersections;
}

static edge_t** insert_edges_and_vtx(vcn_model_t *model,
				     vcn_container_t *vertices,
				     vcn_container_t *edges,
				     bool check_if_edges_exist)
{
	edge_t **edges_array = calloc(model->M, sizeof(*edges_array));
	for (uint32_t i = 0; i < model->M; i++) {
		edge_t aux_edge;
		aux_edge.v1 = insert_vtx_if_not_exist(model, vertices,
						      GET_1_EDGE_VTX(model, i));
		aux_edge.v2 = insert_vtx_if_not_exist(model, vertices,
						      GET_2_EDGE_VTX(model, i));
		edge_set_length(&aux_edge);
		bool allow_to_insert_edge = true;
		if (check_if_edges_exist) {
			edge_t* existing = vcn_container_exist(edges, &aux_edge);
			if (NULL != existing)
				allow_to_insert_edge = false;

		}
		if (allow_to_insert_edge) {
			edge_t *edge = edge_clone(&aux_edge);
			edges_array[i] = edge;
			vcn_container_insert(edges, edge);
		} else {
			edges_array[i] = NULL;
		}
	}
	return edges_array;
}

static vtx_t* insert_vtx_if_not_exist(const vcn_model_t *const model,
				      vcn_container_t *vertices,
				      uint32_t vtx_id)
{
	vtx_t vtx_to_find;
	memcpy(vtx_to_find.x, GET_PVTX(model, vtx_id), 
	       2 * sizeof(*(vtx_to_find.x)));

	vtx_t *vtx = vcn_container_exist(vertices, &vtx_to_find);
	if (NULL == vtx) {
		vtx = vtx_create();
		memcpy(vtx->x, GET_PVTX(model, vtx_id), 2 * sizeof(*(vtx->x)));
		vcn_container_insert(vertices, vtx);
	}
	return vtx;
}

static inline uint32_t intersection_pack_key(const void* const pack)
{
	/* Since status = 3 means a parallel or coincident segments,
	 * the priority to be processed must be high, because it could
	 * dissapear by inserting another intersections. */
  
	void** p = (void**) pack;
	int status = ((int*)p[0])[0];
	return status;
}

static inline edge_t* exist_edge_and_insert_if_not(vcn_container_t* edges, 
						   edge_t* edge)
{
	edge_t* existing_edge = vcn_container_exist(edges, edge);
	if (NULL == existing_edge)
			vcn_container_insert(edges, edge);
	return existing_edge;
}

static void** get_intersection_pack(const edge_t *const sgm1,
				    const edge_t *const sgm2)
{
	/* Calculate intersection */
	double* intersection = calloc(2, sizeof(*intersection));
	int status;
	vcn_utils2D_are_sgm_intersected(sgm1->v1->x,
					sgm1->v2->x,
					sgm2->v1->x,
					sgm2->v2->x,
					intersection, &status);

	/* Verify true intersections */
	if (status == 1 || status == 2) {
		free(intersection);
		return NULL;
	}

	if (status == 3) {
		/* Segments parallel or coincident */
		double area = vcn_utils2D_get_2x_trg_area(sgm1->v1->x,
						 sgm1->v2->x,
						 sgm2->v1->x);
		/* Check if them are parallel */
		if(fabs(area) < NB_GEOMETRIC_TOL){
			/* Collineal segments */
			double length1 = edge_get_length(sgm1);
			double length2 = edge_get_length(sgm2);
			double dist_1A_2A = vcn_utils2D_get_dist(sgm1->v1->x, sgm2->v1->x);
			double dist_1A_2B = vcn_utils2D_get_dist(sgm1->v1->x, sgm2->v2->x);
			double dist_1B_2A = vcn_utils2D_get_dist(sgm1->v2->x, sgm2->v1->x);
			double dist_1B_2B = vcn_utils2D_get_dist(sgm1->v2->x, sgm2->v2->x);
			double max_dist = MAX(dist_1B_2A, dist_1B_2B);
			max_dist = MAX(dist_1A_2B, max_dist);
			max_dist = MAX(dist_1A_2A, max_dist);

			/* Check if they are intersected */
			if(max_dist - (length1 + length2) >= -NB_GEOMETRIC_TOL){
				free(intersection);
				return NULL;
			}
		}else{
			/* Does not consider intersections for parallel segments */
			free(intersection);
			return NULL;
		}
	}
  
	if (status == 4) {
		/* Segments intersecting on sgm1->v1 */
		if(sgm1->v1 == sgm2->v1 || sgm1->v1 == sgm2->v2){
			free(intersection);
			return NULL;
		}
	}

	if (status == 5) {
		/* Segments intersecting on sgm1->v2 */
		if(sgm1->v2 == sgm2->v1 || sgm1->v2 == sgm2->v2){
			free(intersection);
			return NULL;
		}
	}

	if (status == 6) {
		/* Segments intersecting on sgm2->v1 */
		if(sgm2->v1 == sgm1->v1 || sgm2->v1 == sgm1->v2){
			free(intersection);
			return NULL;
		}
	}

	if (status == 7) {
		/* Segments intersecting on sgm2->v2 */
		if(sgm2->v2 == sgm1->v1 || sgm2->v2 == sgm1->v2){
			free(intersection);
			return NULL;
		}
	}

	/* Create intersection pack */
	void **intersection_pack = malloc(4 * sizeof((*intersection_pack)));
	int *intersection_status = malloc(sizeof(*(intersection_status)));
	intersection_status[0] = status;
	intersection_pack[0] = intersection_status;
	intersection_pack[1] = (edge_t*) sgm1;
	intersection_pack[2] = (edge_t*) sgm2;
	intersection_pack[3] = intersection;
  
	/* Return intersection pack */
	return intersection_pack;
}

static void process_segment_intersections(vcn_container_t* avl_sgm,
					  vcn_container_t* sgm_intersect,
					  vcn_container_t* ht_vtx)
{
  

	/* Process intersections */
	while (vcn_container_is_not_empty(sgm_intersect)) {
		void** ipack = (void**)vcn_container_delete_first(sgm_intersect);
		int status = ((int*)ipack[0])[0];
		double* intersection = ipack[3];

		edge_t* sgm1 = ipack[1];
		edge_t* sgm2 = ipack[2];
		if (0 == status) {
			vcn_container_t* sgm_intersect_aux = vcn_container_create(NB_CONTAINER_QUEUE);
			/* Segments intersecting */
			vtx_t* vtx = vtx_create();
			memcpy(vtx->x, intersection, 2 * sizeof(double));
      
			vtx_t* aux_vtx = vcn_container_exist(ht_vtx, vtx);
			if (NULL == aux_vtx) {
				vcn_container_insert(ht_vtx, vtx);
			} else {
				vtx_destroy(vtx);
				vtx = aux_vtx;
			}

			/* Remove segments */
			vcn_container_delete(avl_sgm, sgm1);
			vcn_container_delete(avl_sgm, sgm2);

			/* Insert subsegments of sgm1 */
			edge_t* sgm1B = edge_create();
			sgm1B->v1 = vtx;
			sgm1B->v2 = sgm1->v2;
			sgm1->v2 = vtx;

			/* Update lengths */
			edge_set_length(sgm1);
			edge_set_length(sgm1B);

			/* Insert into the AVL */
			edge_t* existing_sgm1 =
				exist_edge_and_insert_if_not(avl_sgm, sgm1);

			edge_t* existing_sgm1B =
				exist_edge_and_insert_if_not(avl_sgm, sgm1B);
			if (NULL != existing_sgm1B)
				edge_destroy(sgm1B);
			
			void** jpack = search_intersection_pack_with_sgm(sgm_intersect, sgm1);
			while (NULL != jpack) {
				edge_t *sgm_aux = jpack[1];
				if (sgm1 == sgm_aux)
					sgm_aux = jpack[2];
	
				if (NULL == existing_sgm1) {
					void **new_pack = get_intersection_pack(sgm1, sgm_aux);
					if (NULL != new_pack)
						vcn_container_insert(sgm_intersect_aux, new_pack);
				}
	  
				if (NULL == existing_sgm1B) {
					void **new_pack = get_intersection_pack(sgm1B, sgm_aux);
					if (NULL != new_pack)
						vcn_container_insert(sgm_intersect_aux, new_pack);
				}
	
				vcn_container_delete(sgm_intersect, jpack);
	
				destroy_intersection_pack(jpack);
	
				jpack = search_intersection_pack_with_sgm(sgm_intersect, sgm1);
			}

			if(NULL != existing_sgm1)
				edge_destroy(sgm1);

			/* Insert subsegments of sgm2 */
			edge_t* sgm2B = edge_create();
			sgm2B->v1 = vtx;
			sgm2B->v2 = sgm2->v2;
			sgm2->v2 = vtx;

			/* Update lengths */
			edge_set_length(sgm2);
			edge_set_length(sgm2B);

			/* Insert into the AVL */
			edge_t* existing_sgm2 = exist_edge_and_insert_if_not(avl_sgm, sgm2);

			edge_t* existing_sgm2B = exist_edge_and_insert_if_not(avl_sgm, sgm2B);
			if (NULL != existing_sgm2B)
				edge_destroy(sgm2B);
			
			jpack = search_intersection_pack_with_sgm(sgm_intersect, sgm2);
			while (NULL != jpack) {
				edge_t* sgm_aux = jpack[1];
				if (sgm2 == sgm_aux)
					sgm_aux = jpack[2];
	
				if (NULL == existing_sgm2) {
					void** new_pack = get_intersection_pack(sgm2, sgm_aux);
					if (NULL != new_pack)
						vcn_container_insert(sgm_intersect_aux, new_pack);
				}
	
				if (NULL == existing_sgm2B) {
					void** new_pack = get_intersection_pack(sgm2B, sgm_aux);
					if (NULL != new_pack)
						vcn_container_insert(sgm_intersect_aux, new_pack);
				}
	
				vcn_container_delete(sgm_intersect, jpack);
	
				destroy_intersection_pack(jpack);
	
				jpack = search_intersection_pack_with_sgm(sgm_intersect, sgm2);
			}

			if (NULL != existing_sgm2)
				edge_destroy(sgm2);

			while (vcn_container_is_not_empty(sgm_intersect_aux)) {
				void** aux_pack = vcn_container_delete_first(sgm_intersect_aux);
				vcn_container_insert(sgm_intersect, aux_pack);
			}
			vcn_container_destroy(sgm_intersect_aux);
		} else if (3 == status) {
			/* Collineal segments */
			double length1 = edge_get_length(sgm1);
			double length2 = edge_get_length(sgm2);
			double dist_1A_2A = vcn_utils2D_get_dist(sgm1->v1->x, sgm2->v1->x);
			double dist_1A_2B = vcn_utils2D_get_dist(sgm1->v1->x, sgm2->v2->x);
			double dist_1B_2A = vcn_utils2D_get_dist(sgm1->v2->x, sgm2->v1->x);
			double dist_1B_2B = vcn_utils2D_get_dist(sgm1->v2->x, sgm2->v2->x);
			double max_dist = MAX(dist_1B_2A, dist_1B_2B);
			max_dist = MAX(dist_1A_2B, max_dist);
			max_dist = MAX(dist_1A_2A, max_dist);
			/* Verify type of intersection segments */
			if (MAX(length1, length2) - max_dist < -NB_GEOMETRIC_TOL) {
				/* Intersected by one side */
				vcn_container_t* sgm_intersect_aux = vcn_container_create(NB_CONTAINER_QUEUE);
				vcn_container_delete(avl_sgm, sgm1);
				vcn_container_delete(avl_sgm, sgm2);

				if (dist_1A_2A == max_dist)
					sgm1->v2 = sgm2->v1;
				else if (dist_1A_2B == max_dist)
					sgm1->v2 = sgm2->v2;
				else if (dist_1B_2A == max_dist)
					sgm1->v1 = sgm2->v1;
				else
					sgm1->v1 = sgm2->v2;

	    
				/* Update lengths */
				edge_set_length(sgm1);

				/* Insert into the AVL */
				edge_t* existing_sgm1 = exist_edge_and_insert_if_not(avl_sgm, sgm1);

				void** jpack = 
					search_intersection_pack_with_sgm(sgm_intersect, sgm1);
				while (NULL != jpack) {
					edge_t* sgm_aux = (edge_t*)jpack[1];
					if (sgm1 == sgm_aux)
						sgm_aux = (edge_t*)jpack[2];

					if (NULL == existing_sgm1) {
						void** new_pack = get_intersection_pack(sgm1, sgm_aux);
						if (NULL != new_pack)
							vcn_container_insert(sgm_intersect_aux, new_pack);
					}

					vcn_container_delete(sgm_intersect, jpack);
	      
					destroy_intersection_pack(jpack);
	      
					jpack = search_intersection_pack_with_sgm(sgm_intersect, sgm1);
				}
	  
				jpack =
					search_intersection_pack_with_sgm(sgm_intersect, sgm2);
				while (NULL != jpack) {
					edge_t* sgm_aux = jpack[1];
					if(sgm2 == sgm_aux)
						sgm_aux = jpack[2];
					/* Check intersection with segment 1 */
	      
					if(NULL == existing_sgm1) {
						void **kpack = 
							search_intersection_pack_with_both_segments(sgm_intersect_aux, 
												    sgm1, sgm_aux);
						if (NULL == kpack) {
							void **new_pack = get_intersection_pack(sgm1, sgm_aux);
							if (NULL != new_pack)
								vcn_container_insert(sgm_intersect_aux, new_pack);
						}
					}
	      
					vcn_container_delete(sgm_intersect, jpack);
	      
					destroy_intersection_pack(jpack);
	      
					jpack = search_intersection_pack_with_sgm(sgm_intersect, sgm2);
				}	 
	    
				if (NULL != existing_sgm1)
					edge_destroy(sgm1);

				while (vcn_container_is_not_empty(sgm_intersect_aux)) {
					void** aux_pack = vcn_container_delete_first(sgm_intersect_aux);
					vcn_container_insert(sgm_intersect, aux_pack);
				}
				vcn_container_destroy(sgm_intersect_aux);
				edge_destroy(sgm2);
			} else {
				/* One segment contains the other */
				edge_t* small_sgm;
				if (length1 > length2)
					small_sgm = sgm2;  
				else
					small_sgm = sgm1;

				vcn_container_delete(avl_sgm, small_sgm);

				void** jpack =
					search_intersection_pack_with_sgm(sgm_intersect, small_sgm);
				while (NULL != jpack) {
					vcn_container_delete(sgm_intersect, jpack);
					destroy_intersection_pack(jpack);
					jpack = 
						search_intersection_pack_with_sgm(sgm_intersect, small_sgm);
				}
				edge_destroy(small_sgm);
			}
		} else if(4 == status) {
			/* Segments intersecting on sgm1->v1 */
			split_segment_by_vertex(avl_sgm, sgm_intersect, sgm2, sgm1->v1);
    
		} else if(5 == status) {
			/* Segments intersecting on sgm1->v2 */
			split_segment_by_vertex(avl_sgm, sgm_intersect, sgm2, sgm1->v2);
    
		} else if(6 == status) {
			/* Segments intersecting on sgm2->v1 */
			split_segment_by_vertex(avl_sgm, sgm_intersect, sgm1, sgm2->v1);
    
		} else if(7 == status) {
			/* Segments intersecting on sgm2->v2 */
			split_segment_by_vertex(avl_sgm, sgm_intersect, sgm1, sgm2->v2);
		}
		/* Free status */
		destroy_intersection_pack(ipack);
	}
}

static void process_short_edges(vcn_container_t *vertices,
				vcn_container_t *edges,
				double min_length_x_edge)
{
	remove_short_segments(edges, vertices, min_length_x_edge);
	vcn_container_t *intersections = search_intersections_in_edges(edges);
	while (vcn_container_is_not_empty(intersections)) {
		process_segment_intersections(edges, intersections, vertices);
		remove_short_segments(edges, vertices, min_length_x_edge);
		vcn_container_destroy(intersections);
		intersections = search_intersections_in_edges(edges);
	}
	vcn_container_destroy(intersections);
}

static void remove_short_segments(vcn_container_t* avl_sgm,
				  vcn_container_t* ht_vtx,
				  double min_length_x_segment)
{
	while (vcn_container_is_not_empty(avl_sgm)) {
		/* TEMPORAL: This must be faster */
		edge_t* sgm = vcn_container_delete_first(avl_sgm);
		double length = edge_get_length(sgm);
		if (length >= MAX(min_length_x_segment, NB_GEOMETRIC_TOL)) {
			vcn_container_insert(avl_sgm, sgm);
			break;
		}
    
		vcn_container_delete(ht_vtx, sgm->v1);
		vcn_container_delete(ht_vtx, sgm->v2);
		vtx_t* vtx = sgm->v1;

		vcn_container_t* list_sgm = vcn_container_create(NB_CONTAINER_QUEUE);
		vcn_iterator_t* avl_subiter = vcn_iterator_create();
		vcn_iterator_set_container(avl_subiter, avl_sgm);
		while (vcn_iterator_has_more(avl_subiter)) {
			const edge_t* subsgm = vcn_iterator_get_next(avl_subiter);
			if (subsgm->v1 == sgm->v2 || subsgm->v2 == sgm->v2 ||
			    subsgm->v1 == vtx || subsgm->v2 == vtx)
				vcn_container_insert(list_sgm, subsgm);
		}
		vcn_iterator_destroy(avl_subiter);

		vcn_iterator_t* iter = vcn_iterator_create();
		vcn_iterator_set_container(iter, list_sgm);
		while (vcn_iterator_has_more(iter)) {
			edge_t* subsgm = (edge_t*) vcn_iterator_get_next(iter);
			vcn_container_delete(avl_sgm, subsgm);
			if (subsgm->v1 == sgm->v2)
				subsgm->v1 = vtx;
			else if (subsgm->v2 == sgm->v2)
				subsgm->v2 = vtx;
		}
		vcn_iterator_destroy(iter);

		/* Re-insert into the hash tables */
		if (vcn_container_is_not_empty(list_sgm)) {
			vtx->x[0] = (vtx->x[0] + sgm->v2->x[0]) * 0.5;
			vtx->x[1] = (vtx->x[1] + sgm->v2->x[1]) * 0.5;
			vcn_container_insert(ht_vtx, vtx);
		} else {
			vtx_destroy(vtx);
		}

		while (vcn_container_is_not_empty(list_sgm)) {
			edge_t* subsgm = vcn_container_delete_first(list_sgm);
			/* Update length */
			edge_set_length(subsgm);

			/* Insert if it does not exist */
			edge_t* existing_sgm = exist_edge_and_insert_if_not(avl_sgm, subsgm);
			if (NULL != existing_sgm)
				edge_destroy(subsgm);
		}
		vcn_container_destroy(list_sgm);

		/* Free memory */
		vtx_destroy(sgm->v2);
		edge_destroy(sgm);
	}
}

static vcn_container_t* search_intersections_in_edges(vcn_container_t *edges)
{
	vcn_container_t *intersections = 
		vcn_container_create(NB_CONTAINER_SORTED);
	vcn_container_set_key_generator(intersections, intersection_pack_key);
  
	vcn_iterator_t* iter = vcn_iterator_create();
	vcn_iterator_set_container(iter, edges);
	while (vcn_iterator_has_more(iter)) {
		const edge_t* edge1 = vcn_iterator_get_next(iter);
		vcn_iterator_t* subiter = vcn_iterator_clone(iter);
		while (vcn_iterator_has_more(subiter)) {
			const edge_t* edge2 = vcn_iterator_get_next(subiter);
			void** pack = get_intersection_pack(edge1, edge2);
			if (NULL != pack)
				vcn_container_insert(intersections, pack);
		}
		vcn_iterator_destroy(subiter);
	}
	vcn_iterator_destroy(iter);
	return intersections;
}

static void set_as_initial_vtx_in_edges(vcn_container_t *edges)
{
	vcn_iterator_t *iter = vcn_iterator_create();
	vcn_iterator_set_container(iter, edges);
	while (vcn_iterator_has_more(iter)) {
		const edge_t *edge = vcn_iterator_get_next(iter);
		vtx_set_as_initial(edge->v1);
			vtx_set_as_initial(edge->v2);
	}
	vcn_iterator_destroy(iter);
}

static void delete_unused_vertices(vcn_container_t *vertices)
{
	vcn_container_t *to_delete = vcn_container_create(NB_CONTAINER_QUEUE);
	vcn_iterator_t *iter = vcn_iterator_create();
	vcn_iterator_set_container(iter, vertices);
	while (vcn_iterator_has_more(iter)) {
		const vtx_t* vtx = vcn_iterator_get_next(iter);
		if (vtx_is_not_initial(vtx))
			vcn_container_insert(to_delete, vtx);
	}
	vcn_iterator_destroy(iter);

	while (vcn_container_is_not_empty(to_delete)) {
		vtx_t *vtx = vcn_container_delete_first(to_delete);
		vcn_container_delete(vertices, vtx);
		vtx_destroy(vtx);		
	}
	vcn_container_destroy(to_delete);
}

static void** search_intersection_pack_with_sgm
				(const vcn_container_t *const packs,
				 const edge_t *const edge)
{
	vcn_iterator_t* iter = vcn_iterator_create();
	vcn_iterator_set_container(iter, packs);
	void** out_pack = NULL;
	while (vcn_iterator_has_more(iter)) {
		void** pack = (void**) vcn_iterator_get_next(iter);
		if (edge == pack[1] || edge == pack[2]) {
			out_pack = pack;
			break;
		}
	}
	vcn_iterator_destroy(iter);
	return out_pack;
}

static inline void destroy_intersection_pack(void **pack)
{
	free(pack[0]);
	free(pack[3]);
	free(pack);
}

static void** search_intersection_pack_with_both_segments
				(const vcn_container_t *const packs,
				 const edge_t *const edge1,
				 const edge_t *const edge2)
{
	vcn_iterator_t* iter = vcn_iterator_create();
	vcn_iterator_set_container(iter, packs);
	void **out_pack = NULL;
	while (vcn_iterator_has_more(iter)) {
		void** pack = (void**) vcn_iterator_get_next(iter);
		if ((edge1 == pack[1] && edge2 == pack[2]) ||
		    (edge2 == pack[1] && edge1 == pack[2])) {
			out_pack = pack;
			break;
		}
	}
	vcn_iterator_destroy(iter);
	return out_pack;
}

static void split_segment_by_vertex(vcn_container_t* avl_sgm,
				    vcn_container_t* sgm_intersect,
				    edge_t* sgm, vtx_t* vtx)
{
	vcn_container_delete(avl_sgm, sgm);
	edge_t* new_sgm = edge_create();
	
	new_sgm->v1 = vtx;
	new_sgm->v2 = sgm->v2;
	sgm->v2 = vtx;

	/* Update lengths */
	edge_set_length(sgm);
	edge_set_length(new_sgm);

	/* Insert into the AVL */
	edge_t* existing_sgm = exist_edge_and_insert_if_not(avl_sgm, sgm);

	edge_t* existing_new_sgm = exist_edge_and_insert_if_not(avl_sgm, new_sgm);
	if (NULL != existing_new_sgm)
		edge_destroy(new_sgm);
	
	vcn_container_t* sgm_intersect_aux = vcn_container_create(NB_CONTAINER_QUEUE);
  
	void** jpack = search_intersection_pack_with_sgm(sgm_intersect, sgm);
	while (NULL != jpack) {
		edge_t* sgm_aux = jpack[1];
		if (sgm == sgm_aux)
			sgm_aux = jpack[2];

		if (NULL == existing_sgm) {
			void** new_pack = get_intersection_pack(sgm, sgm_aux);
			if (NULL != new_pack)
				vcn_container_insert(sgm_intersect_aux, new_pack);
		}

		if (NULL == existing_new_sgm) {
			void** new_pack = get_intersection_pack(new_sgm, sgm_aux);
			if (NULL != new_pack)
				vcn_container_insert(sgm_intersect_aux, new_pack);
		}
		vcn_container_delete(sgm_intersect, jpack);
		destroy_intersection_pack(jpack);
		jpack = search_intersection_pack_with_sgm(sgm_intersect, sgm);
	}	 

	while (vcn_container_is_not_empty(sgm_intersect_aux)) {
		void** aux_pack = vcn_container_delete_first(sgm_intersect_aux);
		vcn_container_insert(sgm_intersect, aux_pack);
	}
	vcn_container_destroy(sgm_intersect_aux);

	if (NULL != existing_sgm)
		edge_destroy(sgm);
}

vcn_model_t* vcn_model_get_intersection(const vcn_model_t *const model1,
					const vcn_model_t *const model2,
					double min_length_x_segment){
	/* Get combined model */
	vcn_model_t* model = vcn_model_get_combination(model1, model2,
						       min_length_x_segment);
	if(model == NULL) return NULL;

	/* Verify coherence */
	uint32_t model_error_ids[2];
	int model_error = vcn_model_verify_consistence(model, model_error_ids);
	if(model_error != 0){
		/* Report error */
		FILE *fatal_file = fopen("FATAL_ERROR.GEOMETRIC_BOT.log", "a");
		if(fatal_file != NULL){
			char model_file[100];
			void* rand_ptr = malloc(1);
			sprintf(model_file, "FATAL_ERROR.GEOMETRIC_BOT.%p.psl", rand_ptr);
			free(rand_ptr);
			vcn_model_save(model, model_file);
			fprintf(fatal_file, "\nmodel_intersection(...)\n");
			fprintf(fatal_file, "  Incoherent model (%i) -> (%i, %i) \n",
				model_error, model_error_ids[0], model_error_ids[1]);
			fprintf(fatal_file, "  Model saved in %s\n", model_file);
			fclose(fatal_file);
		}
		vcn_model_destroy(model);
		return NULL;
	}

	/* Get holes */
	vcn_mesh_t* mesh = vcn_mesh_create();
	vcn_mesh_generate_from_model(mesh, model);

	uint32_t N_centroids;
	double* centroids = vcn_mesh_get_centroids_of_subareas(mesh, &N_centroids);
	vcn_mesh_destroy(mesh);

	char* mask_centroids = (char*)calloc(N_centroids, 1);

	uint32_t N_new_holes = 0;
	vcn_mesh_t* mesh1 = vcn_mesh_create();
	vcn_mesh_generate_from_model(mesh1, model1);

	vcn_mesh_t* mesh2 = vcn_mesh_create();
	vcn_mesh_generate_from_model(mesh2, model2);
 
	for (uint32_t i=0; i < N_centroids; i++) {
		if (!vcn_mesh_is_vtx_inside(mesh1, &(centroids[i * 2])) ||
		    !vcn_mesh_is_vtx_inside(mesh2, &(centroids[i * 2]))) {
			mask_centroids[i] = 1;
			N_new_holes += 1;
		}
	}
	vcn_mesh_destroy(mesh1);
	vcn_mesh_destroy(mesh2);

	double* new_holes = NULL;
	if (model->H + N_new_holes > 0)
		new_holes =
			(double*) calloc(2 * (model->H + N_new_holes), sizeof(double));
	if (model->H > 0)
		memcpy(new_holes, model->holes, 2 * model->H * sizeof(double));

	N_new_holes = model->H;
	for (uint32_t i=0; i < N_centroids; i++) {
		if (mask_centroids[i] == 0)
			continue;
		memcpy(&(new_holes[N_new_holes * 2]), 
		       &(centroids[i*2]), 2 * sizeof(double));
		N_new_holes += 1;
	}
	free(centroids);
	free(mask_centroids);

	if (model->H > 0)
		free(model->holes);

	model->H = N_new_holes;
	model->holes = new_holes;

	/* Remove isolated segments and vertices */
	mesh = vcn_mesh_create();
	vcn_mesh_generate_from_model(mesh, model);

	vcn_mesh_delete_isolated_segments(mesh);
	vcn_mesh_delete_internal_input_segments(mesh);
	vcn_mesh_delete_isolated_vertices(mesh);

	if (vcn_mesh_get_N_trg(mesh) == 0) {
		vcn_mesh_destroy(mesh);
		vcn_model_destroy(model);
		return NULL;
	}

	vcn_msh3trg_t* msh3trg =
		vcn_mesh_get_msh3trg(mesh, false, true, false, true, true);
	vcn_mesh_destroy(mesh);
  
	vcn_model_destroy(model);
	model = vcn_model_create_from_msh3trg(msh3trg);
	vcn_msh3trg_destroy(msh3trg);

	/* Remove vertices isolated vertices in the interior of the model */
	char* mask = (char*) calloc(model->N, 1);
	for(uint32_t i = 0; i < model->M; i++){
		mask[model->edge[i * 2]] = 1;
		mask[model->edge[i*2+1]] = 1;
	}

	uint32_t N_vtx = 0;
	uint32_t* perm = (uint32_t*) malloc(model->N * sizeof(uint32_t));
	for(uint32_t i = 0; i < model->N; i++){
		if(mask[i] == 0) continue;
		perm[i] = N_vtx;
		N_vtx ++;
	}
	if(N_vtx < model->N){
		double* vertices = (double*) calloc(N_vtx * 2, sizeof(double));
		for(uint32_t i = 0; i < model->N; i++){
			if(mask[i] == 0) continue;
			memcpy(&(vertices[perm[i] * 2]),
			       &(model->vertex[i * 2]), 
			       2 * sizeof(double));
		}  
    
		for(uint32_t i = 0; i < model->M; i++){
			model->edge[i * 2] = perm[model->edge[i * 2]];
			model->edge[i*2+1] = perm[model->edge[i*2+1]];
		}
		model->N = N_vtx;
		free(model->vertex);
		model->vertex = vertices;
	}

	free(perm);
	free(mask);

	/* Return models intersection */
	return model;
}

vcn_model_t* vcn_model_get_union(const vcn_model_t *const model1,
				 const vcn_model_t *const model2,
				 double min_length_x_segment){
  
	/* Get combined model */
	vcn_model_t* model = vcn_model_get_combination(model1, model2,
						       min_length_x_segment);

	if(model == NULL) return NULL;
  
	/* Verify coherence */
	uint32_t model_error_ids[2];
	int model_error = vcn_model_verify_consistence(model, model_error_ids);
	if(model_error != 0){
		/* Report error */
		FILE *fatal_file = fopen("FATAL_ERROR.GEOMETRIC_BOT.log", "a");
		if(fatal_file != NULL){
			char model_file[100];
			void* rand_ptr = malloc(1);
			sprintf(model_file, "FATAL_ERROR.GEOMETRIC_BOT.%p.psl", rand_ptr);
			free(rand_ptr);
			vcn_model_save(model, model_file);
			fprintf(fatal_file, "\nmodel_union(...)\n");
			fprintf(fatal_file, "  Incoherent model (%i) -> (%i, %i) \n",
				model_error, model_error_ids[0], model_error_ids[1]);
			fprintf(fatal_file, "  Model saved in %s\n", model_file);
			fclose(fatal_file);
		}
		vcn_model_destroy(model);
		return NULL;
	}

	/* Remove internal segments and vertices */
	vcn_mesh_t* mesh = vcn_mesh_create();
	vcn_mesh_generate_from_model(mesh, model);

	vcn_mesh_delete_internal_input_segments(mesh);

	vcn_msh3trg_t* msh3trg = vcn_mesh_get_msh3trg(mesh, false, true, false, true, true);
	vcn_mesh_destroy(mesh);
  
	vcn_model_destroy(model);
	model = vcn_model_create_from_msh3trg(msh3trg);
	vcn_msh3trg_destroy(msh3trg);

	/* Remove isolated vertices in the interior of the model */
	char* mask = (char*) calloc(model->N, 1);
	for(uint32_t i = 0; i < model->M; i++){
		mask[model->edge[i * 2]] = 1;
		mask[model->edge[i*2+1]] = 1;
	}

	uint32_t N_vtx = 0;
	uint32_t* perm = (uint32_t*) malloc(model->N * sizeof(uint32_t));
	for(uint32_t i = 0; i < model->N; i++){
		if(mask[i] == 0) continue;
		perm[i] = N_vtx;
		N_vtx ++;
	}
	if(N_vtx < model->N){
		double* vertices = (double*) calloc(N_vtx * 2, sizeof(double));
		for(uint32_t i = 0; i < model->N; i++){
			if(mask[i] == 0) continue;
			memcpy(&(vertices[perm[i] * 2]),
			       &(model->vertex[i * 2]), 
			       2 * sizeof(double));
		}  
    
		for(uint32_t i = 0; i < model->M; i++){
			model->edge[i * 2] = perm[model->edge[i * 2]];
			model->edge[i*2+1] = perm[model->edge[i*2+1]];
		}
		model->N = N_vtx;
		free(model->vertex);
		model->vertex = vertices;
	}

	free(perm);
	free(mask);

	/* Return models intersection */
	return model;
}
