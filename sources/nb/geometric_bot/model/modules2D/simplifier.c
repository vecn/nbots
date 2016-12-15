#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include "nb/memory_bot.h"
#include "nb/container_bot.h"
#include "nb/geometric_bot/utils2D.h"

#include "nb/geometric_bot/model/model2D_struct.h"
#include "nb/geometric_bot/model/model2D.h"
#include "nb/geometric_bot/model/modules2D/simplifier.h"

#define GET_PVTX(model, i) (&((model)->vertex[(i)*2]))
#define GET_1_EDGE_VTX(model, i) ((model)->edge[(i) * 2])
#define GET_2_EDGE_VTX(model, i) ((model)->edge[(i)*2+1])
#define MAX(a, b) (((a)>(b))?(a):(b))
#define POW2(a) ((a)*(a))

static void build_wire(const nb_model_t *const model,
		       nb_container_t *wires, char *mask_wired,
		       uint32_t edge_id);
static bool link_to_wire(const nb_model_t *const model,
			 char *mask_wired, nb_container_t *wire,
			 uint32_t edge_id, uint32_t *vstart,
			 uint32_t *vfinish);
static bool is_at_start(const nb_model_t *const model,
			uint32_t edge_id, uint32_t vstart);
static bool is_at_finish(const nb_model_t *const model,
			 uint32_t edge_id, uint32_t vfinish);
static void link_at_start(const nb_model_t *const model,
			  uint32_t *edge,
			  uint32_t edge_id, uint32_t vstart);
static void insert_at_start(nb_container_t *wire, const void *const val);
static uint32_t update_end_point(const nb_model_t *const model,
				 uint32_t edge_id, uint32_t end_point);
static void link_at_finish(const nb_model_t *const model,
			   uint32_t *edge,
			   uint32_t edge_id, uint32_t vfinish);
static void remove_all_collinearities(const nb_model_t *const model,
				      const nb_container_t *const wires,
				      uint32_t N_fixed_vertices,
				      uint32_t* fixed_vertices,
				      double tolerance);
static void remove_collinearities_of_wire(const nb_model_t *const model,
					  nb_container_t *wire,
					  uint32_t N_fixed_vertices,
					  uint32_t* fixed_vertices,
					  double tolerance);
static bool collapse_wire(const nb_model_t *const model, nb_container_t* wire,
			  uint32_t N_fixed_vertices,
			  uint32_t* fixed_vertices,
			  double tolerance);
static bool collapse_pair_of_edges(const nb_model_t *const model,
				   uint32_t *edge, uint32_t *edge_next,
				   uint32_t N_fixed_vertices,
				   uint32_t* fixed_vertices,
				   double tolerance);

nb_container_t* nb_model_generate_wires(const nb_model_t *const model)
{
	nb_container_t* wires = nb_container_create(NB_QUEUE);
	char* mask_wired = nb_allocate_zero_mem(model->M);
	for (uint32_t i = 0; i < model->M; i++) {
		if (0 == mask_wired[i])
			build_wire(model, wires, mask_wired, i);
	}
	nb_free_mem(mask_wired);
	return wires;
}

static void build_wire(const nb_model_t *const model,
		       nb_container_t *wires, char *mask_wired,
		       uint32_t edge_id)
{
	mask_wired[edge_id] = 1;
	nb_container_t* new_wire = nb_container_create(NB_QUEUE);
	uint32_t* edge = nb_allocate_mem(2 * sizeof(*edge));
	edge[0] = GET_1_EDGE_VTX(model, edge_id);
	edge[1] = GET_2_EDGE_VTX(model, edge_id);
	nb_container_insert(new_wire, edge);
	uint32_t vstart = edge[0];
	uint32_t vfinish = edge[1];
	uint32_t j = 0;
	while (j < model->M) {
		bool linked = false;
		if (0 == mask_wired[j])
			linked = link_to_wire(model, mask_wired, new_wire, j,
					      &vstart, &vfinish);
		if (linked)
			j = 0;
		else
			j += 1;
	}
	nb_container_insert(wires, new_wire);
}

static bool link_to_wire(const nb_model_t *const model,
			 char *mask_wired, nb_container_t *wire,
			 uint32_t edge_id, uint32_t *vstart,
			 uint32_t *vfinish)
{
	bool linked = false;
	if (is_at_start(model, edge_id, *vstart)) {
		mask_wired[edge_id] = 1;
		uint32_t *edge = nb_allocate_mem(2 * sizeof(*edge));
		link_at_start(model, edge, edge_id, *vstart);
		insert_at_start(wire, edge);
		*vstart = update_end_point(model, edge_id, *vstart);
		linked = true;
	} else if (is_at_finish(model, edge_id, *vstart)) {
		mask_wired[edge_id] = 1;
		uint32_t *edge = nb_allocate_mem(2 * sizeof(*edge));
		link_at_finish(model, edge, edge_id, *vfinish);
		nb_container_insert(wire, edge);
		*vfinish = update_end_point(model, edge_id, *vfinish);
		linked = true;
	}
	return linked;
}

static inline bool is_at_start(const nb_model_t *const model,
			       uint32_t edge_id, uint32_t vstart)
{
	return vstart == GET_1_EDGE_VTX(model, edge_id) ||
		vstart == GET_2_EDGE_VTX(model, edge_id);
}

static inline bool is_at_finish(const nb_model_t *const model,
				uint32_t edge_id, uint32_t vfinish)
{
	return vfinish == GET_1_EDGE_VTX(model, edge_id) ||
		vfinish == GET_2_EDGE_VTX(model, edge_id);
}

static inline void link_at_start(const nb_model_t *const model,
				 uint32_t *edge,
				 uint32_t edge_id, uint32_t vstart)
{	
	if (vstart == GET_1_EDGE_VTX(model, edge_id)) {
		edge[0] = GET_2_EDGE_VTX(model, edge_id);
		edge[1] = GET_1_EDGE_VTX(model, edge_id);
	} else {
		edge[0] = GET_1_EDGE_VTX(model, edge_id);
		edge[1] = GET_2_EDGE_VTX(model, edge_id);
	}
}

static inline void insert_at_start(nb_container_t *wire,
				   const void *const val)
{
	int8_t status;
	nb_container_do(wire, "insert_first", (void*) val, &status);
	assert(0 == status);
}

static inline uint32_t update_end_point(const nb_model_t *const model,
					uint32_t edge_id, uint32_t end_point)
{
	uint32_t new_end_point;
	if (end_point == GET_1_EDGE_VTX(model, edge_id))
		new_end_point = GET_2_EDGE_VTX(model, edge_id);
	else
		new_end_point = GET_1_EDGE_VTX(model, edge_id);
	return new_end_point;

}

static inline void link_at_finish(const nb_model_t *const model,
				  uint32_t *edge,
				  uint32_t edge_id, uint32_t vfinish)
{
	if (vfinish == GET_1_EDGE_VTX(model, edge_id)) {
		edge[0] = GET_1_EDGE_VTX(model, edge_id);
		edge[1] = GET_2_EDGE_VTX(model, edge_id);
	} else {
		edge[0] = GET_2_EDGE_VTX(model, edge_id);
		edge[1] = GET_1_EDGE_VTX(model, edge_id);
	}
}

void nb_model_collapse_small_segments(nb_model_t* model,
				       double tolerance,
				       uint32_t N_fixed_vertices,
				       uint32_t* fixed_vertices)
/* WARNING: Could produce segment intersections */
{
	double tol2 = POW2(tolerance);
	while (true) {
		/* Mark fixed vertices */
		char* mask_fixed_vertices = nb_allocate_zero_mem(model->N);
		for (uint32_t i = 0; i < N_fixed_vertices; i++)
			mask_fixed_vertices[fixed_vertices[i]] = 1;
    
		/* Allocate masks for removed vertices */
		char* mask_vtx_removed = nb_allocate_zero_mem(model->N);

		/* Remove segments */
		uint32_t N_removed = 0;
		for (uint32_t i = 0; i < model->M; /* Nothing here */) {
			uint32_t v1 = model->edge[i * 2];
			uint32_t v2 = model->edge[i*2+1];
			if (mask_fixed_vertices[v1] && mask_fixed_vertices[v2]) {
				i++;
				continue;
			}
    
			double dist2 = nb_utils2D_get_dist2(&(model->vertex[v1*2]),
						    &(model->vertex[v2*2]));
			if (dist2 < tol2) {
				N_removed += 1;
				/* Remove segment */
				if (mask_fixed_vertices[v2]) {
					uint32_t aux = v1;
					v1 = v2;
					v2 = aux;
				} else if (!mask_fixed_vertices[v1]) {
					model->vertex[v1 * 2] = 0.5 *
						(model->vertex[v1 * 2] + model->vertex[v2 * 2]);
					model->vertex[v1*2+1] = 0.5 *
						(model->vertex[v1*2+1] + model->vertex[v2*2+1]);
				}
				mask_vtx_removed[v2] = 1;
				/* (TEMPORAL: Could be faster using tables) */
				for (uint32_t j = 0; j < model->M; j++) {
					if (i == j)
						continue;
					if (model->edge[j * 2] == v2)
						model->edge[j * 2] = v1;
					if (model->edge[j*2+1] == v2)
						model->edge[j*2+1] = v1;
				}

				/* Remove from array */
				uint32_t N_remaining = model->M - i - 1;
				if (N_remaining > 0) {
					uint32_t* remaining_sgm =
						nb_allocate_mem(N_remaining * 2 * sizeof(uint32_t));
					memcpy(remaining_sgm, &(model->edge[(i+1)*2]),
					       N_remaining * 2 * sizeof(uint32_t));
					memcpy(&(model->edge[i*2]), remaining_sgm,
					       N_remaining * 2 * sizeof(uint32_t));
					nb_free_mem(remaining_sgm);
				}
				model->M -= 1;
			} else {
				i++;
			}
		}

		if (N_removed == 0) 
			break;

		/* Remove repeated segments (TEMPORAL: Could be faster) */
		for (uint32_t i = 0; i < model->M; /* Nothing here */) {
			uint32_t v1 = model->edge[i * 2];
			uint32_t v2 = model->edge[i*2+1];
			bool is_repeated = false;
			for (uint32_t j = 0; j < model->M; j++) {
				if (i == j)
					continue;
				if ((model->edge[j * 2] == v1 && model->edge[j*2+1] == v2) ||
				    (model->edge[j * 2] == v2 && model->edge[j*2+1] == v1)){
					is_repeated = true;
					break;
				}
			}
			if (is_repeated) {
				/* Remove from array */
				uint32_t N_remaining = model->M - i - 1;
				if (N_remaining > 0) {
					uint32_t* remaining_sgm =
						nb_allocate_mem(N_remaining * 2 * sizeof(uint32_t));
					memcpy(remaining_sgm, &(model->edge[(i+1)*2]),
					       N_remaining * 2 * sizeof(uint32_t));
					memcpy(&(model->edge[i*2]), remaining_sgm,
					       N_remaining * 2 * sizeof(uint32_t));
					nb_free_mem(remaining_sgm);
				}
				model->M -= 1;
			} else {
				i++;
			}
		}

		/* Allocate new arrays */
		uint32_t N_vtx = model->N - N_removed;
		double* vertices = nb_allocate_mem(N_vtx * 2 * sizeof(*vertices));
		uint32_t* segments = nb_allocate_mem(model->M * 2 * sizeof(*segments));

		memcpy(segments, model->edge, model->M * 2 * sizeof(uint32_t));
		nb_free_mem(model->edge);
		model->edge = segments;

		uint32_t* perm_vtx = nb_allocate_zero_mem(model->N *
							  sizeof(*perm_vtx));
		uint32_t vtx_counter = 0;
		for (uint32_t i = 0; i < model->N; i++) {
			if (mask_vtx_removed[i])
				continue;
			memcpy(&(vertices[vtx_counter*2]),
			       &(model->vertex[i*2]), 
			       2 * sizeof(double));
			perm_vtx[i] = vtx_counter;
			vtx_counter += 1;
		}
		nb_free_mem(model->vertex);
		model->vertex = vertices;
		model->N = N_vtx;

		for (uint32_t i = 0; i < model->M; i++) {
			model->edge[i * 2] = perm_vtx[model->edge[i * 2]];
			model->edge[i*2+1] = perm_vtx[model->edge[i*2+1]];
		}
		for (uint32_t i = 0; i < N_fixed_vertices; i++)
			fixed_vertices[i] = perm_vtx[fixed_vertices[i]];
		nb_free_mem(perm_vtx);

		/* Free memory */
		nb_free_mem(mask_fixed_vertices);
		nb_free_mem(mask_vtx_removed);
	}
}


void nb_model_collapse_colinear_vertices(nb_model_t* model,
					  uint32_t N_fixed_vertices,
					  uint32_t* fixed_vertices,
					  double tolerance){
	nb_container_t* wires = nb_model_generate_wires(model);

	remove_all_collinearities(model, wires, N_fixed_vertices,
				  fixed_vertices, tolerance);

	/* Count surviving vertices and segments */
	uint32_t N_sgm = 0;
	uint32_t N_vtx = 0;
	char* mask_vtx = nb_allocate_zero_mem(model->N);
	uint32_t* perm_vtx = nb_allocate_zero_mem(model->N * sizeof(uint32_t));
	nb_iterator_t *iter = nb_iterator_create();
	nb_iterator_set_container(iter, wires);
	while (nb_iterator_has_more(iter)) {
		const nb_container_t* wire = nb_iterator_get_next(iter);
		N_sgm += nb_container_get_length(wire);    
		nb_iterator_t* subiter = nb_iterator_create();
		nb_iterator_set_container(subiter, wire);
		while (nb_iterator_has_more(subiter)) {
			const uint32_t* sgm = nb_iterator_get_next(subiter);
			if (!mask_vtx[sgm[0]]) {
				mask_vtx[sgm[0]] = 1;
				perm_vtx[sgm[0]] = N_vtx++;
			}
			if (!mask_vtx[sgm[1]]) {
				mask_vtx[sgm[1]] = 1;
				perm_vtx[sgm[1]] = N_vtx++;
			}
		}
		nb_iterator_destroy(subiter);
	}
	nb_iterator_destroy(iter);

	/* Set new vertices and edges */
	double* vertices = nb_allocate_mem(2 * N_vtx * sizeof(*vertices));
	uint32_t* segments = nb_allocate_mem(2 * N_sgm * sizeof(*segments));

	uint32_t isgm = 0;
	while (nb_container_is_not_empty(wires)) {
		nb_container_t* wire = nb_container_delete_first(wires);
		while (nb_container_is_not_empty(wire)) {
			uint32_t* sgm = nb_container_delete_first(wire);
			segments[isgm * 2] = perm_vtx[sgm[0]];
			segments[isgm*2+1] = perm_vtx[sgm[1]];
			isgm++;
			memcpy(&(vertices[perm_vtx[sgm[0]] * 2]),
			       &(model->vertex[sgm[0] * 2]),
			       2 * sizeof(double));
			memcpy(&(vertices[perm_vtx[sgm[1]] * 2]),
			       &(model->vertex[sgm[1] * 2]),
			       2 * sizeof(double));
		}
		nb_container_destroy(wire);
	}
	nb_free_mem(perm_vtx);
	nb_free_mem(mask_vtx);

	model->N = N_vtx;
	nb_free_mem(model->vertex);
	model->vertex = vertices;

	model->M = N_sgm;
	nb_free_mem(model->edge);
	model->edge = segments;

	/* Free memory */
	nb_container_destroy(wires);
}

static void remove_all_collinearities(const nb_model_t *const model,
				      const nb_container_t *const wires,
				      uint32_t N_fixed_vertices,
				      uint32_t* fixed_vertices,
				      double tolerance)
{
	nb_iterator_t* iter = nb_iterator_create();
	nb_iterator_set_container(iter, wires);
	while (nb_iterator_has_more(iter)) {
		const nb_container_t* wire = nb_iterator_get_next(iter);
		remove_collinearities_of_wire(model, (nb_container_t*)wire,
					      N_fixed_vertices, fixed_vertices,
					      tolerance);
	}
	nb_iterator_destroy(iter);
}

static inline void remove_collinearities_of_wire(const nb_model_t *const model,
						 nb_container_t *wire,
						 uint32_t N_fixed_vertices,
						 uint32_t* fixed_vertices,
						 double tolerance)
{
	bool collapsed = true;
	while (collapsed)
	  collapsed = collapse_wire(model, wire, N_fixed_vertices,
				    fixed_vertices, tolerance);
}

static bool collapse_wire(const nb_model_t *const model,
			  nb_container_t* wire,
			  uint32_t N_fixed_vertices,
			  uint32_t* fixed_vertices,
			  double tolerance)
{
	nb_container_t* clone = nb_container_clone(wire);
	nb_container_clear(wire);
	uint32_t* edge = nb_container_delete_first(clone);
	bool collapsed = false;
	while (nb_container_is_not_empty(clone)) {
		uint32_t* edge_next = nb_container_delete_first(clone);
		collapsed = collapse_pair_of_edges(model, edge, edge_next,
						   N_fixed_vertices,
						   fixed_vertices,
						   tolerance);
		if (!collapsed) {
			nb_container_insert(wire, edge);
			edge = edge_next;
		}
	}
	nb_container_insert(wire, edge);
	nb_container_destroy(clone);
	return collapsed;
}

static bool collapse_pair_of_edges(const nb_model_t *const model,
				   uint32_t *edge, uint32_t *edge_next,
				   uint32_t N_fixed_vertices,
				   uint32_t* fixed_vertices,
				   double tolerance)
{
	double tol2 = POW2(tolerance);
	uint32_t v1;
	uint32_t v_middle;
	uint32_t v2;
	if (edge[0] == edge_next[0]) {
		v1 = edge[1];
		v_middle = edge[0];
		v2 = edge_next[1];
	} else if (edge[0] == edge_next[1]) {
		v1 = edge[1];
		v_middle = edge[0];
		v2 = edge_next[0];
	} else if (edge[1] == edge_next[0]) {
		v1 = edge[0];
		v_middle = edge[1];
		v2 = edge_next[1];
	} else {
		v1 = edge[0];
		v_middle = edge[1];
		v2 = edge_next[0];
	}
	bool will_be_removed = false;
	if (nb_utils2D_pnt_lies_on_sgm(&(model->vertex[v1*2]),
						 &(model->vertex[v2*2]),
						 &(model->vertex[v_middle*2])))
		will_be_removed = true;
	
	bool is_fixed = false;
	for (uint32_t i = 0; i < N_fixed_vertices; i++) {
		if (fixed_vertices[i] == v_middle) {
			is_fixed = true;
			break;
		}
	}

	if (!will_be_removed && !is_fixed) {
		double closest_point[2];
		nb_utils2D_get_closest_pnt_to_sgm(&(model->vertex[v1*2]),
						   &(model->vertex[v2*2]),
						   &(model->vertex[v_middle*2]),
						   closest_point);
		double dist_p2 = nb_utils2D_get_dist2(&(model->vertex[v_middle*2]),
						       closest_point);
		if (dist_p2 < tol2)
			will_be_removed = true;
	}

	bool collapsed = false;
	if (will_be_removed && !is_fixed) {
		nb_free_mem(edge_next);
		edge[0] = v1;
		edge[1] = v2;
		collapsed = true;
	}
	return collapsed;
}

void nb_model_unify_edge(nb_model_t* model, double* vtx1, double* vtx2)
/* If the edge is not unified yet, it would have the last index */
{
	uint32_t v1 = model->N;
	uint32_t v2 = model->N;
	uint32_t N_vtx_to_remove = 0;
	char* mask_vtx_to_remove = nb_allocate_zero_mem(model->N);
	for (uint32_t i = 0; i < model->N; i++) {
		if (v1 == model->N) {
			if (nb_utils2D_get_dist2(vtx1, &(model->vertex[i*2])) < NB_GEOMETRIC_TOL_POW2) {
				v1 = i;
				continue;
			}
		}
		if (v2 == model->N) {
			if (nb_utils2D_get_dist2(vtx2, &(model->vertex[i*2])) < NB_GEOMETRIC_TOL_POW2) {
				v2 = i;
				continue;
			}
		}
  
		if (nb_utils2D_pnt_lies_on_sgm(vtx1, vtx2,
						&(model->vertex[i * 2]))) {
			N_vtx_to_remove += 1;
			mask_vtx_to_remove[i] = 1;
		}
	}

	uint32_t N_sgm_to_remove = 0;
	char* mask_sgm_to_remove = nb_allocate_zero_mem(model->M);
	for (uint32_t i = 0; i < model->M; i++) {
		uint32_t sv1 = model->edge[i * 2];
		uint32_t sv2 = model->edge[i*2+1];
		if((sv1 == v1 && sv2 == v2) ||
		   (sv1 == v2 && sv2 == v1)){
			nb_free_mem(mask_vtx_to_remove);
			nb_free_mem(mask_sgm_to_remove);
			return;
		}
		if (mask_vtx_to_remove[sv1] ||
		   mask_vtx_to_remove[sv2]) {
			N_sgm_to_remove += 1;
			mask_sgm_to_remove[i] = 1;
		}
	}
	/* Allocate new arrays */
	double* vertices = 
		nb_allocate_mem((model->N - N_vtx_to_remove) * 2 * sizeof(*vertices));
	uint32_t* edges =
		nb_allocate_mem((model->M - N_sgm_to_remove + 1) * 2 * sizeof(*edges));

	/* Fill new array with vertices */
	uint32_t vtx_counter = 0;
	uint32_t* perm_vtx = nb_allocate_mem(model->N * sizeof(*perm_vtx));
	for (uint32_t i = 0; i < model->N; i++) {
		if (mask_vtx_to_remove[i])
			continue;
		memcpy(&(vertices[vtx_counter * 2]),
		       &(model->vertex[i * 2]),
		       2 * sizeof(double));
		perm_vtx[i] = vtx_counter;
		vtx_counter += 1;
	}

	/* Fill new array with edges */
	uint32_t sgm_counter = 0;
	for (uint32_t i = 0; i < model->M; i++) {
		if (mask_sgm_to_remove[i])
			continue;
		edges[sgm_counter * 2] = perm_vtx[model->edge[i * 2]];
		edges[sgm_counter*2+1] = perm_vtx[model->edge[i*2+1]];
		sgm_counter += 1;
	}
	edges[sgm_counter * 2] = perm_vtx[v1];
	edges[sgm_counter*2+1] = perm_vtx[v2];

	/* Update data in the model */
	model->N = model->N - N_vtx_to_remove;
	model->M = model->M - N_sgm_to_remove + 1;
	nb_free_mem(model->vertex);
	model->vertex = vertices;
	nb_free_mem(model->edge);
	model->edge = edges;

	/* Free memory */
	nb_free_mem(perm_vtx);
	nb_free_mem(mask_vtx_to_remove);
	nb_free_mem(mask_sgm_to_remove);
}
