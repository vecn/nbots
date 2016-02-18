#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>
#include <math.h>

#include "nb/math_bot.h"
#include "nb/container_bot/container.h"
#include "nb/container_bot/iterator.h"
#include "nb/geometric_bot/utils2D.h"
#include "nb/geometric_bot/mesh/mesh2D.h"

#include "../model2D_struct.h"
#include "nb/geometric_bot/model/model2D.h"
#include "nb/geometric_bot/model/modules2D/verifier.h"

#define GET_PVTX(model, i) (&((model)->vertex[(i)*2]))
#define GET_1_EDGE_VTX(model, i) ((model)->edge[(i) * 2])
#define GET_2_EDGE_VTX(model, i) ((model)->edge[(i)*2+1])
#define MAX(a, b) (((a)>(b))?(a):(b))

static bool vtx_is_out_of_bounds(const vcn_model_t *const model,
				 uint32_t vtx_id,
				 uint32_t ids_of_error[2]);

int vcn_model_verify_consistence(const vcn_model_t *const model,
				 uint32_t ids_causing_error[2])
{
	int error_code;
	if (!vcn_model_have_vertices(model))
		error_code = 1;
	else if (!vcn_model_have_edges(model))
		error_code = 2;
	else if (vcn_model_have_repeated_vertices(model, ids_causing_error))
		error_code = 3;
	else if (vcn_model_have_incoherent_edges(model, ids_causing_error))
		error_code = 4;
	else if (vcn_model_have_repeated_edges(model, ids_causing_error))
		error_code = 5;
	else if (vcn_model_have_intersected_edges(model, ids_causing_error))
		error_code = 6;
	else if (vcn_model_have_vtx_intersecting_edges(model, ids_causing_error))
		error_code = 7;
	else if (vcn_model_have_unclosed_boundary(model))
		error_code = 8;
	else
		error_code = 0;

	return error_code;
}


inline bool vcn_model_have_vertices(const vcn_model_t *const model)
{
	return (0 < model->N);
}

inline bool vcn_model_have_edges(const vcn_model_t *const model)
{
	return (0 < model->M);
}

bool vcn_model_have_repeated_vertices(const vcn_model_t *const model,
				      uint32_t repeated_ids[2])
{
	bool have_rep_vtx = false;
	for (uint32_t i = 0; i < model->N && !have_rep_vtx; i++) {
		for(uint32_t j = i+1; j < model->N; j++) {
			double dist = vcn_utils2D_get_dist2(GET_PVTX(model, i),
						   GET_PVTX(model, j));
			if (dist < NB_GEOMETRIC_TOL) {
				if (NULL != repeated_ids) {
					repeated_ids[0] = i;
					repeated_ids[1] = j;
				}
				have_rep_vtx = true;
				break;
			}	
		}
	}
	return have_rep_vtx;
}

bool vcn_model_have_incoherent_edges(const vcn_model_t *const model,
				     uint32_t ids_edge_and_vtx[2])
{
	bool have_inc_edge = false;
	for (uint32_t i=0; i < model->M; i++) {
		bool v1_out = vtx_is_out_of_bounds(model, 
						   GET_1_EDGE_VTX(model, i),
						   ids_edge_and_vtx);
		bool v2_out = vtx_is_out_of_bounds(model, 
						   GET_2_EDGE_VTX(model, i),
						   ids_edge_and_vtx);
		if (v1_out || v2_out) {
			if (NULL != ids_edge_and_vtx)
				ids_edge_and_vtx[0] = i;
			have_inc_edge = true;
			break;
		}
	}
	return have_inc_edge;
}

static bool vtx_is_out_of_bounds(const vcn_model_t *const model,
				 uint32_t vtx_id,
				 uint32_t ids_of_error[2])
{
	bool out = false;
	if (vtx_id >= model->N) {
		if (NULL != ids_of_error)
			ids_of_error[1] = vtx_id;
		out = true;
	}
	return out;
}

bool vcn_model_have_repeated_edges(const vcn_model_t *const model,
				   uint32_t repeated_ids[2])
{
	bool have_rep_sgm = false;
	for (uint32_t i = 0; i < model->M && !have_rep_sgm; i++) {
		for (uint32_t j = i+1; j < model->M; j++) {
			bool equals_i1_and_j1 = 
				GET_1_EDGE_VTX(model, i) == GET_1_EDGE_VTX(model, j);
			bool equals_i2_and_j2 = 
				GET_2_EDGE_VTX(model, i) == GET_2_EDGE_VTX(model, j);
			bool equals_i1_and_j2 = 
				GET_1_EDGE_VTX(model, i) == GET_2_EDGE_VTX(model, j);
			bool equals_i2_and_j1 = 
				GET_2_EDGE_VTX(model, i) == GET_1_EDGE_VTX(model, j);
			bool flag_A = equals_i1_and_j1 && equals_i2_and_j2;
			bool flag_B = equals_i1_and_j2 && equals_i2_and_j1;
			if (flag_A || flag_B) {
				if (NULL != repeated_ids) {
					repeated_ids[0] = i;
					repeated_ids[1] = j;
				}
				have_rep_sgm = true;
				break;
			}
		}
	}
	return have_rep_sgm;
}

bool vcn_model_have_intersected_edges(const vcn_model_t *const model,
				      uint32_t intersected_ids[2])
{
	bool out = false;
	for (uint32_t i = 0; i < model->M-1 && !out; i++) {
		for (uint32_t j = i+1; j < model->M; j++) {
			uint32_t v1 = GET_1_EDGE_VTX(model, i);
			uint32_t v2 = GET_2_EDGE_VTX(model, i);
			uint32_t v3 = GET_1_EDGE_VTX(model, j);
			uint32_t v4 = GET_2_EDGE_VTX(model, j);

			nb_intersect_t status = 
				vcn_utils2D_are_sgm_intersected
						(GET_PVTX(model, v1),
						 GET_PVTX(model, v2),
						 GET_PVTX(model, v3),
						 GET_PVTX(model, v4),
						 NULL);
			bool are_intersecting = NB_INTERSECTED == status;
      
			if (NB_PARALLEL == status ||
			    NB_INTERSECT_ON_A1 == status ||
			    NB_INTERSECT_ON_A2 == status ||
			    NB_INTERSECT_ON_B1 == status ||
			    NB_INTERSECT_ON_B2 == status) {
				double dist_1A_2A = vcn_utils2D_get_dist(GET_PVTX(model, v1),
							     GET_PVTX(model, v3));
				double dist_1A_2B = vcn_utils2D_get_dist(GET_PVTX(model, v1),
							     GET_PVTX(model, v4));
				double dist_1B_2A = vcn_utils2D_get_dist(GET_PVTX(model, v2),
							     GET_PVTX(model, v3));
				double dist_1B_2B = vcn_utils2D_get_dist(GET_PVTX(model, v2),
							     GET_PVTX(model, v4));

				if (NB_PARALLEL == status) {
					/* Segments parallel or coincident */
					double area = vcn_utils2D_get_2x_trg_area(GET_PVTX(model, v1),
									 GET_PVTX(model, v2),
									 GET_PVTX(model, v3));
					/* Check if them are parallel */
					if (fabs(area) < NB_GEOMETRIC_TOL) {
						/* Collineal segments */
						double length1 = vcn_utils2D_get_dist(GET_PVTX(model, v1),
									  GET_PVTX(model, v2));
						double length2 = vcn_utils2D_get_dist(GET_PVTX(model, v3),
									  GET_PVTX(model, v4));
						double max_dist = MAX(dist_1B_2A, dist_1B_2B);
						max_dist = MAX(max_dist, dist_1A_2B);
						max_dist = MAX(max_dist, dist_1A_2A);
						/* Check if they are intersected */
						if(max_dist - (length1 + length2) < 
						   -NB_GEOMETRIC_TOL)
							are_intersecting = true;
					}
				}
				if (NB_INTERSECT_ON_A1 == status) {
					/* Segments intersecting on sgm1->v1 */
					if(dist_1A_2A > NB_GEOMETRIC_TOL &&
					   dist_1A_2B > NB_GEOMETRIC_TOL)
						are_intersecting = true;
				}

				if (NB_INTERSECT_ON_A2 == status) {
					/* Segments intersecting on sgm1->v2 */
					if(dist_1B_2A > NB_GEOMETRIC_TOL &&
					   dist_1B_2B > NB_GEOMETRIC_TOL)
						are_intersecting = true;
				}

				if (NB_INTERSECT_ON_B1 == status) {
					/* Segments intersecting on sgm2->v1 */
					if(dist_1A_2A > NB_GEOMETRIC_TOL &&
					   dist_1B_2A > NB_GEOMETRIC_TOL)
						are_intersecting = true;
				}

				if (NB_INTERSECT_ON_B2 == status) {
					/* Segments intersecting on sgm2->v2 */
					if(dist_1A_2B > NB_GEOMETRIC_TOL &&
					   dist_1B_2B > NB_GEOMETRIC_TOL)
						are_intersecting = true;
				}
			}
       
			if (are_intersecting) {
				if (NULL != intersected_ids) {
					intersected_ids[0] = i;
					intersected_ids[1] = j;
				}
				out = true;
				break;
			}
		}
	}
	return out;
}

bool vcn_model_have_vtx_intersecting_edges(const vcn_model_t *const model,
					   uint32_t ids_edge_and_vtx[2])
{
	bool out = false;
	for (uint32_t i = 0; i < model->M && !out; i++) {
		for (uint32_t j = 0; j < model->N; j++) {
			if (model->edge[i * 2] == j)
				continue;
			if (model->edge[i*2+1] == j)
				continue;
			uint32_t id1 = model->edge[i * 2];
			uint32_t id2 = model->edge[i*2+1];
			double x1 = model->vertex[id1 * 2];
			double y1 = model->vertex[id1*2+1];
			double x2 = model->vertex[id2 * 2];
			double y2 = model->vertex[id2*2+1];
			double x3 = model->vertex[j * 2];
			double y3 = model->vertex[j*2+1];
			double a1 = (x2-x3)/(x2-x1);
			double a2 = (y2-y3)/(y2-y1);
			if (fabs(a1-a2) < NB_GEOMETRIC_TOL) {
				if (a1 > 0 && a1 < 1) {
					if (NULL != ids_edge_and_vtx) {
						ids_edge_and_vtx[0] = i; /* Edge id */
						ids_edge_and_vtx[1] = j; /* Vertex id */
					}
					out = true;
					break;
				}
			}
		}
	}
	return out;
}

bool vcn_model_have_unclosed_boundary(const vcn_model_t *const model)
{
	/* Verify unclosed shapes and unknown errors*/
	vcn_mesh_t* mesh = vcn_mesh_create();
	vcn_mesh_generate_from_model(mesh, model);

	bool is_unclosed;
	if (vcn_mesh_is_empty(mesh))
		is_unclosed = true;

	vcn_mesh_destroy(mesh);

	return is_unclosed;
}
