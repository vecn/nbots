#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>
#include <math.h>
#include <alloca.h>

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

static void* ipack_create(void);
static void ipack_destroy(void *ipack_ptr);
static int8_t ipack_comparer(const void *ipack1_ptr,
			     const void *ipack2_ptr);

static double get_scale_and_disp(const vcn_model_t *const model1,
				 const vcn_model_t *const model2,
				 double *xdisp, double *ydisp);
static void scale_and_displace(vcn_model_t *model, double scale,
			       double xdisp, double ydisp);
static nb_container_t* search_intersections(vcn_model_t *model1,
					     vcn_model_t *model2,
					     nb_container_t *vertices,
					     nb_container_t *edges);
static edge_t** insert_edges_and_vtx(vcn_model_t *model,
				      nb_container_t *vertices,
				      nb_container_t *edges,
				      bool check_if_edges_exist);
static vtx_t* insert_vtx_if_not_exist(const vcn_model_t *const model,
				      nb_container_t *vertices,
				      uint32_t vtx_id);
static edge_t* exist_edge_and_insert_if_not(nb_container_t* edges,
					    edge_t* edge);
static ipack_t* get_intersection_pack(const edge_t *const sgm1,
				      const edge_t *const sgm2);
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
static void calculate_collinear_dists(const ipack_t *ipack,
				      collinear_dists *cdists);
static bool collinear_are_not_contained(collinear_dists *cdists);
static void process_collinear_intersected(ipack_t *ipack,
					  collinear_dists *cdists,
					  nb_container_t* avl_sgm,
					  nb_container_t* sgm_intersect);
static void process_collinear_contained(ipack_t *ipack,
					collinear_dists *cdists,
					nb_container_t* avl_sgm,
					nb_container_t* sgm_intersect);
static void process_short_edges(nb_container_t *vertices,
				nb_container_t *edges,
				double min_length_x_edg);
static void remove_short_segments(nb_container_t* avl_sgm,
				  nb_container_t* ht_vtx,
				  double min_length_x_segment);
static nb_container_t* search_intersections_in_edges(nb_container_t *edges);
static void set_as_initial_vtx_in_edges(nb_container_t *edges);
static void delete_unused_vertices(nb_container_t *vertices);
static ipack_t* search_intersection_pack_with_sgm
				(const nb_container_t *const packs,
				 const edge_t *const edge);
static ipack_t* search_intersection_pack_with_both_segments
			(const nb_container_t *const packs,
			 const edge_t *const edge1,
			 const edge_t *const edge2);
static void split_segment_by_vertex(nb_container_t* avl_sgm,
				    nb_container_t* sgm_intersect,
				    edge_t* sgm, vtx_t* vtx);
static void set_intersection_from_combination(vcn_model_t *model,
					      const vcn_model_t *const model1,
					      const vcn_model_t *const model2);
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
static void set_union_from_combination(vcn_model_t *model);
static void set_substraction_from_combination(vcn_model_t *model,
					      const vcn_model_t *const model2);
static void set_substraction_holes(vcn_model_t *model,
				   const vcn_model_t *const model2);
static uint32_t mask_substraction_holes(const vcn_model_t *model2,
					uint32_t N_centroids,
					const double *centroids,
					char *mask_centroids);

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
	nb_container_t* ht_vtx = nb_container_create(NB_HASH);
	nb_container_set_key_generator(ht_vtx, vtx_hash_key);
	nb_container_set_comparer(ht_vtx, vtx_compare);
	nb_container_set_destroyer(ht_vtx, vtx_destroy);

	nb_container_t* avl_sgm = nb_container_create(NB_SORTED);
	nb_container_set_comparer(avl_sgm, edge_compare);
	nb_container_set_destroyer(avl_sgm, edge_destroy);
	   
	nb_container_t* sgm_intersect = search_intersections(model1, model2,
							     ht_vtx, avl_sgm);

	process_segment_intersections(avl_sgm, sgm_intersect, ht_vtx);

	nb_container_destroy(sgm_intersect);

	process_short_edges(ht_vtx, avl_sgm, min_length_x_segment * scale);

	/* Check that a volume is possible */
	if (nb_container_get_length(avl_sgm) < 3) {
		vcn_model_destroy(model1);
		vcn_model_destroy(model2);
		nb_container_destroy(avl_sgm);
		nb_container_destroy(ht_vtx);
		return vcn_model_create();
	}
	
	set_as_initial_vtx_in_edges(avl_sgm);
	delete_unused_vertices(ht_vtx);

	/* Allocate new model */
	vcn_model_t* model = vcn_model_create();
	model->N = nb_container_get_length(ht_vtx);
	model_alloc_vertices(model);
	model->M = nb_container_get_length(avl_sgm);
	model_alloc_edges(model);

	/* Set vertices to the model */
	uint32_t vtx_id = 0;
	nb_iterator_t *ht_iter = nb_iterator_create();
	nb_iterator_set_container(ht_iter, ht_vtx);
	while(nb_iterator_has_more(ht_iter)){
		const vtx_t* vtx = nb_iterator_get_next(ht_iter);
		vtx_set_id((vtx_t*)vtx, vtx_id);
		memcpy(GET_PVTX(model, vtx_id), vtx->x, 2 * sizeof(*(vtx->x)));
		vtx_id++;
	}
	nb_iterator_destroy(ht_iter);

	/* Set segments to the model */
	uint32_t sgm_counter = 0;
	nb_iterator_t *avl_iter = nb_iterator_create();
	nb_iterator_set_container(avl_iter, avl_sgm);
	while(nb_iterator_has_more(avl_iter)){
		const edge_t* edge = nb_iterator_get_next(avl_iter);
		model->edge[sgm_counter * 2] = vtx_get_id(edge->v1);
		model->edge[sgm_counter*2+1] = vtx_get_id(edge->v2);
		sgm_counter++;
	}
	nb_iterator_destroy(avl_iter);

	/* Get holes */
	vcn_mesh_t* mesh = vcn_mesh_create();
	vcn_mesh_set_geometric_constraint(mesh, NB_MESH_GEOM_CONSTRAINT_MIN_ANGLE, 0);
	vcn_mesh_generate_from_model(mesh, model);

	uint32_t N_centroids;
	double* centroids = vcn_mesh_get_centroids_of_subareas(mesh, &N_centroids);
	vcn_mesh_destroy(mesh);

	if (N_centroids > 1) {
		char* mask_centroids = calloc(N_centroids, 1);
		
		vcn_mesh_t* mesh1 = vcn_mesh_create();
		vcn_mesh_set_geometric_constraint(mesh1, NB_MESH_GEOM_CONSTRAINT_MIN_ANGLE, 0);
		vcn_mesh_generate_from_model(mesh1, model1);

		vcn_mesh_t* mesh2 = vcn_mesh_create();
		vcn_mesh_set_geometric_constraint(mesh2, NB_MESH_GEOM_CONSTRAINT_MIN_ANGLE, 0);
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
			for (uint32_t i = 0; i < N_centroids; i++) {
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
  
	nb_container_set_destroyer(ht_vtx, vtx_destroy);
	nb_container_destroy(ht_vtx);
  
	while (nb_container_is_not_empty(avl_sgm)) {
		edge_t* sgm = nb_container_delete_first(avl_sgm);
		edge_destroy(sgm);
	}
	nb_container_destroy(avl_sgm);

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

static void* ipack_create(void)
{
	return calloc(1, sizeof(ipack_t));
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

static nb_container_t* search_intersections(vcn_model_t *model1,
					     vcn_model_t *model2,
					     nb_container_t *vertices,
					     nb_container_t *edges)
{
	edge_t** edges1 = insert_edges_and_vtx(model1, vertices, edges, false);
	edge_t** edges2 = insert_edges_and_vtx(model2, vertices, edges, true);
	nb_container_t *intersections =	nb_container_create(NB_SORTED);
	nb_container_set_comparer(intersections, ipack_comparer);
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
	free(edges1);
	free(edges2);
	return intersections;
}

static edge_t** insert_edges_and_vtx(vcn_model_t *model,
				     nb_container_t *vertices,
				     nb_container_t *edges,
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
	return edges_array;
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
	ipack_t *ipack = NULL;
	/* Calculate intersection */
	double intersection[2];
	nb_intersect_t status =
		vcn_utils2D_are_sgm_intersected(sgm1->v1->x,
						sgm1->v2->x,
						sgm2->v1->x,
						sgm2->v2->x,
						intersection);

	/* Verify true intersections */
	if (NB_NOT_INTERSECTED == status)
		goto EXIT;

	if (NB_PARALLEL == status) {
		/* Segments parallel or coincident */
		double area = vcn_utils2D_get_2x_trg_area(sgm1->v1->x,
						 sgm1->v2->x,
						 sgm2->v1->x);
		/* Check if them are parallel */
		if (fabs(area) < NB_GEOMETRIC_TOL) {
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
			if (max_dist - (length1 + length2) >= -NB_GEOMETRIC_TOL)
				goto EXIT;
		} else {
			goto EXIT;
		}
	}
  
	if (NB_INTERSECT_ON_A1 == status) {
		/* Segments intersecting on sgm1->v1 */
		if (sgm1->v1 == sgm2->v1 || sgm1->v1 == sgm2->v2)
			goto EXIT;
	}

	if (NB_INTERSECT_ON_A2 == status) {
		/* Segments intersecting on sgm1->v2 */
		if (sgm1->v2 == sgm2->v1 || sgm1->v2 == sgm2->v2)
			goto EXIT;
	}

	if (NB_INTERSECT_ON_B1 == status) {
		/* Segments intersecting on sgm2->v1 */
		if (sgm2->v1 == sgm1->v1 || sgm2->v1 == sgm1->v2)
			goto EXIT;
	}

	if (NB_INTERSECT_ON_B2 == status) {
		/* Segments intersecting on sgm2->v2 */
		if (sgm2->v2 == sgm1->v1 || sgm2->v2 == sgm1->v2)
			goto EXIT;
	}

	/* Create intersection pack */
	ipack = ipack_create();
	ipack->status = status;
	ipack->sgm1 = (edge_t*) sgm1;
	ipack->sgm2 = (edge_t*) sgm2;
	memcpy(ipack->intersection, intersection, 2 * sizeof(*intersection));
EXIT:
	return ipack;
}

static void process_segment_intersections(nb_container_t* avl_sgm,
					  nb_container_t* sgm_intersect,
					  nb_container_t* ht_vtx)
{
	/* Process intersections */
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
	case NB_NOT_INTERSECTED:
		/* Do nothing */
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
	calculate_collinear_dists(ipack, &cdists);
	if (collinear_are_not_contained(&cdists))
		process_collinear_intersected(ipack, &cdists,
					      avl_sgm, sgm_intersect);
	else
		process_collinear_contained(ipack, &cdists, avl_sgm,
					    sgm_intersect);
}

static void calculate_collinear_dists(const ipack_t *ipack,
				      collinear_dists *cdists)
{
	cdists->length1 = edge_get_length(ipack->sgm1);
	cdists->length2 = edge_get_length(ipack->sgm2);
	cdists->dist_1A_2A = vcn_utils2D_get_dist(ipack->sgm1->v1->x,
						  ipack->sgm2->v1->x);
	cdists->dist_1A_2B = vcn_utils2D_get_dist(ipack->sgm1->v1->x,
						  ipack->sgm2->v2->x);
	cdists->dist_1B_2A = vcn_utils2D_get_dist(ipack->sgm1->v2->x,
						  ipack->sgm2->v1->x);
	cdists->dist_1B_2B = vcn_utils2D_get_dist(ipack->sgm1->v2->x,
						  ipack->sgm2->v2->x);
	cdists->max_dist = MAX(cdists->dist_1B_2A, cdists->dist_1B_2B);
	cdists->max_dist = MAX(cdists->dist_1A_2B, cdists->max_dist);
	cdists->max_dist = MAX(cdists->dist_1A_2A, cdists->max_dist);
}

static inline bool collinear_are_not_contained(collinear_dists *cdists)
{
	double max_length = MAX(cdists->length1, cdists->length2);
	return (max_length - cdists->max_dist < -NB_GEOMETRIC_TOL);

}

static void process_collinear_intersected(ipack_t *ipack,
					  collinear_dists *cdists,
					  nb_container_t* avl_sgm,
					  nb_container_t* sgm_intersect)
{
	/* Intersected by one side */
	nb_container_t* sgm_intersect_aux = nb_container_create(NB_QUEUE);
	nb_container_delete(avl_sgm, ipack->sgm1);
	nb_container_delete(avl_sgm, ipack->sgm2);

	if (cdists->dist_1A_2A == cdists->max_dist)
		ipack->sgm1->v2 = ipack->sgm2->v1;
	else if (cdists->dist_1A_2B == cdists->max_dist)
		ipack->sgm1->v2 = ipack->sgm2->v2;
	else if (cdists->dist_1B_2A == cdists->max_dist)
		ipack->sgm1->v1 = ipack->sgm2->v1;
	else
		ipack->sgm1->v1 = ipack->sgm2->v2;

	    
	/* Update lengths */
	edge_set_length(ipack->sgm1);

	/* Insert into the AVL */
	edge_t* existing_sgm1 =
		exist_edge_and_insert_if_not(avl_sgm, ipack->sgm1);

	ipack_t* jpack = 
		search_intersection_pack_with_sgm(sgm_intersect,
						  ipack->sgm1);
	while (NULL != jpack) {
		edge_t* sgm_aux = jpack->sgm1;
		if (ipack->sgm1 == sgm_aux)
			sgm_aux = jpack->sgm2;

		if (NULL == existing_sgm1) {
			ipack_t* new_pack =
				get_intersection_pack(ipack->sgm1,
						      sgm_aux);
			if (NULL != new_pack)
				nb_container_insert(sgm_intersect_aux,
						    new_pack);
		}

		nb_container_delete(sgm_intersect, jpack);
	      
		ipack_destroy(jpack);
	      
		jpack = search_intersection_pack_with_sgm(sgm_intersect,
							  ipack->sgm1);
	}
	jpack = search_intersection_pack_with_sgm(sgm_intersect,
						  ipack->sgm2);
	while (NULL != jpack) {
		edge_t* sgm_aux = jpack->sgm1;
		if (ipack->sgm2 == sgm_aux)
			sgm_aux = jpack->sgm2;
		/* Check intersection with segment 1 */
	      
		if(NULL == existing_sgm1) {
			ipack_t *kpack = 
				search_intersection_pack_with_both_segments(sgm_intersect_aux,
									    ipack->sgm1, sgm_aux);
			if (NULL == kpack) {
				ipack_t *new_pack =
					get_intersection_pack(ipack->sgm1,
							      sgm_aux);
				if (NULL != new_pack)
					nb_container_insert(sgm_intersect_aux, new_pack);
			}
		}
	      
		nb_container_delete(sgm_intersect, jpack);
	      
		ipack_destroy(jpack);
	      
		jpack = search_intersection_pack_with_sgm(sgm_intersect, ipack->sgm2);
	}	 
	    
	if (NULL != existing_sgm1)
		edge_destroy(ipack->sgm1);

	while (nb_container_is_not_empty(sgm_intersect_aux)) {
		void** aux_pack = nb_container_delete_first(sgm_intersect_aux);
		nb_container_insert(sgm_intersect, aux_pack);
	}
	nb_container_destroy(sgm_intersect_aux);
	edge_destroy(ipack->sgm2);
}

static void process_collinear_contained(ipack_t *ipack,
					collinear_dists *cdists,
					nb_container_t* avl_sgm,
					nb_container_t* sgm_intersect)
{
	edge_t* small_sgm;
	if (cdists->length1 > cdists->length2)
		small_sgm = ipack->sgm2;  
	else
		small_sgm = ipack->sgm1;

	nb_container_delete(avl_sgm, small_sgm);

	ipack_t* jpack =
		search_intersection_pack_with_sgm(sgm_intersect, 
						  small_sgm);
	while (NULL != jpack) {
		nb_container_delete(sgm_intersect, jpack);
		ipack_destroy(jpack);
		jpack = search_intersection_pack_with_sgm(sgm_intersect,
							  small_sgm);
	}
	edge_destroy(small_sgm);
}

static void process_short_edges(nb_container_t *vertices,
				nb_container_t *edges,
				double min_length_x_edge)
{
	remove_short_segments(edges, vertices, min_length_x_edge);
	nb_container_t *intersections = search_intersections_in_edges(edges);
	while (nb_container_is_not_empty(intersections)) {
		process_segment_intersections(edges, intersections, vertices);
		remove_short_segments(edges, vertices, min_length_x_edge);
		nb_container_destroy(intersections);
		intersections = search_intersections_in_edges(edges);
	}
	nb_container_destroy(intersections);
}

static void remove_short_segments(nb_container_t* avl_sgm,
				  nb_container_t* ht_vtx,
				  double min_length_x_segment)
{
	while (nb_container_is_not_empty(avl_sgm)) {
		/* TEMPORAL: This must be faster */
		edge_t* sgm = nb_container_delete_first(avl_sgm);
		double length = edge_get_length(sgm);
		if (length >= MAX(min_length_x_segment, NB_GEOMETRIC_TOL)) {
			nb_container_insert(avl_sgm, sgm);
			break;
		}
    
		nb_container_delete(ht_vtx, sgm->v1);
		nb_container_delete(ht_vtx, sgm->v2);
		vtx_t* vtx = sgm->v1;

		nb_container_t* list_sgm = nb_container_create(NB_QUEUE);
		nb_iterator_t* avl_subiter = nb_iterator_create();
		nb_iterator_set_container(avl_subiter, avl_sgm);
		while (nb_iterator_has_more(avl_subiter)) {
			const edge_t* subsgm = nb_iterator_get_next(avl_subiter);
			if (subsgm->v1 == sgm->v2 || subsgm->v2 == sgm->v2 ||
			    subsgm->v1 == vtx || subsgm->v2 == vtx)
				nb_container_insert(list_sgm, subsgm);
		}
		nb_iterator_destroy(avl_subiter);

		nb_iterator_t* iter = nb_iterator_create();
		nb_iterator_set_container(iter, list_sgm);
		while (nb_iterator_has_more(iter)) {
			edge_t* subsgm = (edge_t*) nb_iterator_get_next(iter);
			nb_container_delete(avl_sgm, subsgm);
			if (subsgm->v1 == sgm->v2)
				subsgm->v1 = vtx;
			else if (subsgm->v2 == sgm->v2)
				subsgm->v2 = vtx;
		}
		nb_iterator_destroy(iter);

		/* Re-insert into the hash tables */
		if (nb_container_is_not_empty(list_sgm)) {
			vtx->x[0] = (vtx->x[0] + sgm->v2->x[0]) * 0.5;
			vtx->x[1] = (vtx->x[1] + sgm->v2->x[1]) * 0.5;
			nb_container_insert(ht_vtx, vtx);
		} else {
			vtx_destroy(vtx);
		}

		while (nb_container_is_not_empty(list_sgm)) {
			edge_t* subsgm = nb_container_delete_first(list_sgm);
			/* Update length */
			edge_set_length(subsgm);

			/* Insert if it does not exist */
			edge_t* existing_sgm = exist_edge_and_insert_if_not(avl_sgm, subsgm);
			if (NULL != existing_sgm)
				edge_destroy(subsgm);
		}
		nb_container_destroy(list_sgm);

		/* Free memory */
		vtx_destroy(sgm->v2);
		edge_destroy(sgm);
	}
}

static nb_container_t* search_intersections_in_edges(nb_container_t *edges)
{
	nb_container_t *intersections = nb_container_create(NB_SORTED);
	nb_container_set_comparer(intersections, ipack_comparer);
  
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
	return intersections;
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

static ipack_t* search_intersection_pack_with_sgm
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

static ipack_t* search_intersection_pack_with_both_segments
				(const nb_container_t *const packs,
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
	
	nb_container_t* sgm_intersect_aux = nb_container_create(NB_QUEUE);
  
	ipack_t* ipack = search_intersection_pack_with_sgm(sgm_intersect, sgm);
	while (NULL != ipack) {
		edge_t* sgm_aux = ipack->sgm1;
		if (sgm == sgm_aux)
			sgm_aux = ipack->sgm2;

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

		ipack = search_intersection_pack_with_sgm(sgm_intersect, sgm);
	}	 

	while (nb_container_is_not_empty(sgm_intersect_aux)) {
		ipack_t* aux_pack = nb_container_delete_first(sgm_intersect_aux);
		nb_container_insert(sgm_intersect, aux_pack);
	}
	nb_container_destroy(sgm_intersect_aux);

	if (NULL != existing_sgm)
		edge_destroy(sgm);
}

vcn_model_t* vcn_model_get_intersection(const vcn_model_t *const model1,
					const vcn_model_t *const model2,
					double min_length_x_segment)
{
	vcn_model_t* model = vcn_model_get_combination(model1, model2,
						       min_length_x_segment);
	if (NULL != model)
		set_intersection_from_combination(model, model1, model2);
	return model;

}

static void set_intersection_from_combination(vcn_model_t *model,
					      const vcn_model_t *const model1,
					      const vcn_model_t *const model2)
{
	set_intersection_holes(model, model1, model2);
	delete_isolated_elements(model);
	delete_isolated_internal_vtx(model);
}

static void set_intersection_holes(vcn_model_t *model,
				   const vcn_model_t *const model1,
				   const vcn_model_t *const model2)
{
	uint32_t N_centroids;
	double *centroids = 
		get_centroids_of_model_subareas(model, &N_centroids);

	char* mask_centroids = alloca(N_centroids);

	uint32_t N_new_holes =
		mask_intersection_holes(model1, model2, N_centroids,
					centroids, mask_centroids);
	
	set_new_holes_to_model(N_new_holes, N_centroids,
			       centroids, mask_centroids, model);
	free(centroids);
}

static double* get_centroids_of_model_subareas(const vcn_model_t *model,
						uint32_t *N_centroids)
{
	vcn_mesh_t* mesh = vcn_mesh_create();
	vcn_mesh_set_geometric_constraint(mesh, NB_MESH_GEOM_CONSTRAINT_MIN_ANGLE, 0);
	vcn_mesh_generate_from_model(mesh, model);
	double* centroids = 
		vcn_mesh_get_centroids_of_subareas(mesh, N_centroids);
	vcn_mesh_destroy(mesh);
	return centroids;
}

static uint32_t mask_intersection_holes(const vcn_model_t *model1,
					const vcn_model_t *model2,
					uint32_t N_centroids,
					const double *centroids,
					char *mask_centroids)
{
	uint32_t N_new_holes = 0;
	vcn_mesh_t* mesh1 = vcn_mesh_create();
	vcn_mesh_set_geometric_constraint(mesh1, NB_MESH_GEOM_CONSTRAINT_MIN_ANGLE, 0);
	vcn_mesh_generate_from_model(mesh1, model1);

	vcn_mesh_t* mesh2 = vcn_mesh_create();
	vcn_mesh_set_geometric_constraint(mesh2, NB_MESH_GEOM_CONSTRAINT_MIN_ANGLE, 0);
	vcn_mesh_generate_from_model(mesh2, model2);
 
	for (uint32_t i = 0; i < N_centroids; i++) {
		if (!vcn_mesh_is_vtx_inside(mesh1, &(centroids[i * 2])) ||
		    !vcn_mesh_is_vtx_inside(mesh2, &(centroids[i * 2]))) {
			mask_centroids[i] = 1;
			N_new_holes += 1;
		} else {
			mask_centroids[i] = 0;
		}
	}
	vcn_mesh_destroy(mesh1);
	vcn_mesh_destroy(mesh2);
	return N_new_holes;
}

static void set_new_holes_to_model(uint32_t N_new_holes,
				   uint32_t N_centroids,
				   const double *centroids,
				   const char *mask_centroids,
				   vcn_model_t *model)
{
	if (0 < N_new_holes) {
		double* new_holes = calloc(2 * (model->H + N_new_holes),
					   sizeof(*new_holes));
		if (model->H > 0)
			memcpy(new_holes, model->holes,
			       2 * model->H * sizeof(*new_holes));

		N_new_holes = model->H;
		for (uint32_t i = 0; i < N_centroids; i++) {
			if (1 == mask_centroids[i]) {
				memcpy(&(new_holes[N_new_holes * 2]), 
				       &(centroids[i*2]),
				       2 * sizeof(*new_holes));
				N_new_holes += 1;
			}
		}

		if (model->H > 0)
			free(model->holes);

		model->H = model->H + N_new_holes;
		model->holes = new_holes;
	}
}

static void delete_isolated_elements(vcn_model_t *model)
{
	vcn_mesh_t *mesh = vcn_mesh_create();
	vcn_mesh_set_geometric_constraint(mesh, NB_MESH_GEOM_CONSTRAINT_MIN_ANGLE, 0);
	vcn_mesh_generate_from_model(mesh, model);

	vcn_mesh_delete_isolated_segments(mesh);
	vcn_mesh_delete_internal_input_segments(mesh);
	vcn_mesh_delete_isolated_vertices(mesh);

	if (0 == vcn_mesh_get_N_trg(mesh)) {
		vcn_mesh_destroy(mesh);
		vcn_model_clear(model);
		return;
	}

	vcn_msh3trg_t* msh3trg =
		vcn_mesh_get_msh3trg(mesh, false, true, false, true, true);
	vcn_mesh_destroy(mesh);

	vcn_model_generate_from_msh3trg(model, msh3trg);
	vcn_msh3trg_destroy(msh3trg);
}

static void delete_isolated_internal_vtx(vcn_model_t *model)
{
	char* mask = calloc(model->N, 1);
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
}

vcn_model_t* vcn_model_get_union(const vcn_model_t *model1,
				 const vcn_model_t *model2,
				 double min_length_x_segment){
  
	/* Get combined model */
	vcn_model_t* model = vcn_model_get_combination(model1, model2,
						       min_length_x_segment);

	if (NULL != model)
		set_union_from_combination(model);
	return model;
}

static void set_union_from_combination(vcn_model_t *model)
{
	delete_isolated_elements(model);
	delete_isolated_internal_vtx(model);
}

vcn_model_t* vcn_model_get_substraction(const vcn_model_t *const model1,
					const vcn_model_t *const model2,
					double min_length_x_segment)
{
	vcn_model_t* model = vcn_model_get_combination(model1, model2,
						       min_length_x_segment);
	if (NULL != model)
		set_substraction_from_combination(model, model2);
	return model;

}

static void set_substraction_from_combination(vcn_model_t *model,
					      const vcn_model_t *const model2)
{
	set_substraction_holes(model, model2);
	delete_isolated_elements(model);
	delete_isolated_internal_vtx(model);	
}


static void set_substraction_holes(vcn_model_t *model,
				   const vcn_model_t *const model2)
{
	uint32_t N_centroids;
	double *centroids = 
		get_centroids_of_model_subareas(model, &N_centroids);

	char* mask_centroids = alloca(N_centroids);
	
	uint32_t N_new_holes =
		mask_substraction_holes(model2, N_centroids,
					centroids, mask_centroids);

	set_new_holes_to_model(N_new_holes, N_centroids,
			       centroids, mask_centroids, model);

	free(centroids);	
}

static uint32_t mask_substraction_holes(const vcn_model_t *model2,
					uint32_t N_centroids,
					const double *centroids,
					char *mask_centroids)
{
	uint32_t N_new_holes = 0;

	vcn_mesh_t* mesh2 = vcn_mesh_create();
	vcn_mesh_set_geometric_constraint(mesh2, NB_MESH_GEOM_CONSTRAINT_MIN_ANGLE, 0);
	vcn_mesh_generate_from_model(mesh2, model2);
 
	for (uint32_t i = 0; i < N_centroids; i++) {
		if (vcn_mesh_is_vtx_inside(mesh2, &(centroids[i * 2]))) {
			mask_centroids[i] = 1;
			N_new_holes += 1;
		} else {
			mask_centroids[i] = 0;
		}
	}
	vcn_mesh_destroy(mesh2);
	return N_new_holes;
}
