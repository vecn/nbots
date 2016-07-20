#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>
#include <math.h>
#include <alloca.h>

#include "nb/math_bot.h"
#include "nb/memory_bot.h"
#include "nb/container_bot/container.h"
#include "nb/container_bot/iterator.h"
#include "nb/geometric_bot/utils2D.h"
#include "nb/geometric_bot/mesh/mesh2D.h"
#include "nb/geometric_bot/mesh/modules2D/area_analizer.h"

#include "nb/geometric_bot/model/model2D_struct.h"
#include "nb/geometric_bot/model/model2D.h"

#define GET_PVTX(model, i) (&((model)->vertex[(i)*2]))
#define GET_1_EDGE_VTX(model, i) ((model)->edge[(i) * 2])
#define GET_2_EDGE_VTX(model, i) ((model)->edge[(i)*2+1])
#define MAX(a, b) (((a)>(b))?(a):(b))

static inline void *malloc_model(void);

inline uint16_t vcn_model_get_memsize(void)
{
	return sizeof(vcn_model_t);
}

inline void vcn_model_init(void *model_ptr)
{
	memset(model_ptr, 0, vcn_model_get_memsize());
}

void vcn_model_copy(void *model_ptr, const void *src_model_ptr)
{
	vcn_model_t *model = model_ptr;
	const vcn_model_t *src_model = src_model_ptr;

	model->N = src_model->N;
	if (0 < model->N) {
		uint32_t size = 2 * model->N * sizeof(*(model->vertex));
		model->vertex = malloc(size);
		memcpy(model->vertex, src_model->vertex, size);
	}

	model->M = src_model->M;
	if (0 < model->M) {
		uint32_t size = 2 * model->M * sizeof(*(model->edge));
		model->edge = malloc(size);
		memcpy(model->edge, src_model->edge, size);
	}

	model->H = src_model->H;
	if (0 < model->H) {
		uint32_t size = 2 * model->H * sizeof(*(model->holes));
		model->holes = malloc(size);
		memcpy(model->holes, src_model->holes, size);
	}
}

inline void vcn_model_finish(void *model_ptr)
{
	vcn_model_clear(model_ptr);
}

inline void* vcn_model_create(void)
{
	vcn_model_t *model = malloc_model();
	vcn_model_init(model);
	return model;
}

static inline void *malloc_model(void)
{
	uint32_t size = vcn_model_get_memsize();
	return malloc(size);
}

inline void* vcn_model_clone(const void *model_ptr)
{
	vcn_model_t *model = malloc_model();
	vcn_model_copy(model, model_ptr);
	return model;
}

inline void vcn_model_destroy(void *model_ptr)
{
	vcn_model_finish(model_ptr);
	free(model_ptr);
}

void vcn_model_clear(void *model_ptr)
{
	vcn_model_t *model = model_ptr;
	if (0 < model->N)
		free(model->vertex);
	model->N = 0;
	model->vertex = NULL;
	
	if (0 < model->M)
		free(model->edge);
	model->M = 0;
	model->edge = NULL;
	
	if (0 < model->H)
		free(model->holes);
	model->H = 0;
	model->holes = NULL;
}

vcn_model_t* vcn_model_load(const char* filename)
{
	/* OPPORTUNITY: Use custom format */
	vcn_model_t* model = vcn_model_create();
	/* Open file and verify if exist */
	FILE* fp = fopen(filename, "r");
	if (NULL == fp)
		return NULL;
	/* Read number of vertices */
	if (1 != fscanf(fp, "%u", &(model->N)))
		return NULL;
	if (0 == model->N)
		return NULL;
	nb_model_alloc_vertices(model);
	/* Read vertices */
	for (int32_t i = 0; i < model->N; i++) {
		if (fscanf(fp, "%lf %lf",
			   &(model->vertex[i * 2]),
			   &(model->vertex[i*2+1])) != 2)
			return NULL;
	}
	/* Read number of edge */
	if (fscanf(fp, "%u", &(model->M)) != 1)
		return NULL;

	nb_model_alloc_edges(model);
	/* Read segments */
	for (int32_t i = 0; i < model->M; i++) {
	  if (fscanf(fp, "%u %u",
		     &(model->edge[i * 2]),
		     &(model->edge[i*2+1])) != 2)
	    return NULL;
	}

	/* Read number of holes */
	if (fscanf(fp, "%u", &(model->H)) != 1)
		return NULL;

	nb_model_alloc_holes(model);
	for (int32_t i = 0; i < model->H; i++) {
	  if (fscanf(fp, "%lf %lf",
		     &(model->holes[i * 2]),
		     &(model->holes[i*2+1])) != 2)
	    return NULL;
	}
	/* Close file */
	fclose(fp);
	/* Successful operation */
	return model;
}

void vcn_model_load_from_farrays(vcn_model_t *model,
				 float vertices[], uint32_t N_vertices,
				 uint32_t segments[], uint32_t N_segments,
				 float holes[], uint32_t N_holes)
{
	model->N = N_vertices;
	model->M = N_segments;
	model->H = N_holes;

	nb_model_alloc_vertices(model);
	for (uint32_t i = 0; i < 2 * N_vertices; i++)
	  model->vertex[i] = vertices[i];

	nb_model_alloc_edges(model);
	for (uint32_t i = 0; i < 2 * N_segments; i++)
	  model->edge[i] = segments[i];

	nb_model_alloc_holes(model);
	for (uint32_t i = 0; i < 2 * N_holes; i++)
	  model->holes[i] = holes[i];
}

vcn_model_t* vcn_model_create_rectangle(double x_min, double y_min,
					double x_max, double y_max)
{
	vcn_model_t* model = vcn_model_create();
	model->N = 4;
	model->M = 4;
	nb_model_alloc_vertices(model);
	nb_model_alloc_edges(model);
  
	/* Set vertices */
	model->vertex[0] = x_min;
	model->vertex[1] = y_min;
	model->vertex[2] = x_max;
	model->vertex[3] = y_min;
	model->vertex[4] = x_max;
	model->vertex[5] = y_max;
	model->vertex[6] = x_min;
	model->vertex[7] = y_max;

	/* Set edges */
	model->edge[0] = 0;
	model->edge[1] = 1;
	model->edge[2] = 1;
	model->edge[3] = 2;
	model->edge[4] = 2;
	model->edge[5] = 3;
	model->edge[6] = 3;
	model->edge[7] = 0;

	return model;
}

vcn_model_t* vcn_model_create_polygon(double radius,
				      double x_center,
				      double y_center,
				      uint32_t N_sides)
{
	vcn_model_t *model = vcn_model_create();
	double angle_step = (NB_MATH_PI * 2.0) / N_sides;
	model->N = N_sides;
	model->M = N_sides;
	nb_model_alloc_vertices(model);
	nb_model_alloc_edges(model);
	for (uint32_t i = 0; i < model->N; i++) {
		double angle = i * angle_step;
		model->vertex[i * 2] = x_center + radius * cos(angle);
		model->vertex[i*2+1] = y_center + radius * sin(angle);
		model->edge[i * 2] = i;
		model->edge[i*2+1] = (i+1) % model->N;
	}
	return model;
}

vcn_model_t* vcn_model_create_circle(double radius,
				     double x_center,
				     double y_center,
				     double side_length)
{
	double perimeter = 2.0 * NB_MATH_PI * radius;
	uint32_t n = (uint32_t) (perimeter / side_length + 0.5);
	if (10 > n)
		n = 10;
	return vcn_model_create_polygon(radius, x_center, y_center, n);
}

vcn_model_t* vcn_model_create_from_msh3trg
		(const vcn_msh3trg_t *const restrict msh3trg)
{
	vcn_model_t *model = vcn_model_create();
	vcn_model_generate_from_msh3trg(model, msh3trg);
	return model;
}

void vcn_model_generate_from_msh3trg(vcn_model_t *model,
				     const vcn_msh3trg_t *const msh3trg)
{
	uint32_t sgm_memsize = 2 * msh3trg->N_input_segments * sizeof(uint32_t);
	uint32_t* segments = NB_SOFT_MALLOC(sgm_memsize);

	uint32_t vtx_memsize = 2 * msh3trg->N_input_vertices * sizeof(double);
	double* vertices = NB_SOFT_MALLOC(vtx_memsize);

	uint32_t idx_memsize = msh3trg->N_input_vertices * sizeof(uint32_t);
	uint32_t* vtx_index_relation = NB_SOFT_MALLOC(idx_memsize);

	uint32_t N_vertices = 0;
	for (uint32_t i = 0; i < msh3trg->N_input_vertices; i++) {
		if (msh3trg->input_vertices[i] < msh3trg->N_vertices) {
			uint32_t idx = msh3trg->input_vertices[i];
			vertices[N_vertices * 2] = msh3trg->vertices[idx * 2];
			vertices[N_vertices*2+1] = msh3trg->vertices[idx*2+1];
			vtx_index_relation[i] = N_vertices;
			N_vertices += 1;
		}
	}

	uint32_t N_segments = 0;
	for (uint32_t i = 0; i < msh3trg->N_input_segments; i++) {
		if (0 < msh3trg->N_subsgm_x_inputsgm[i]) {
			uint32_t v1 = msh3trg->meshvtx_x_inputsgm[i][0];
			uint32_t last_idx = msh3trg->N_subsgm_x_inputsgm[i];
			uint32_t v2 = msh3trg->meshvtx_x_inputsgm[i][last_idx];
			for (uint32_t j = 0; j < msh3trg->N_input_vertices; j++) {
				if (msh3trg->input_vertices[j] == v1) {
					v1 = j;
					break;
				}
			}
			for (uint32_t j = 0; j < msh3trg->N_input_vertices; j++) {
				if (msh3trg->input_vertices[j] == v2) {
					v2 = j;
					break;
				}
			}
			segments[N_segments * 2] = vtx_index_relation[v1];
			segments[N_segments*2+1] = vtx_index_relation[v2];
			N_segments += 1;
		}
	}
	/* Build model without holes */
	vcn_model_clear(model);
	model->N = N_vertices;
	model->vertex = malloc(2 * N_vertices * sizeof(*(model->vertex)));
	memcpy(model->vertex, vertices, 
	       2 * N_vertices * sizeof(*(model->vertex)));
	model->M = N_segments;
	model->edge = malloc(2 * N_segments * sizeof(*(model->edge)));
	memcpy(model->edge, segments, 2 * N_segments * sizeof(*(model->edge)));

	NB_SOFT_FREE(sgm_memsize, segments);
	NB_SOFT_FREE(vtx_memsize, vertices);
	NB_SOFT_FREE(idx_memsize, vtx_index_relation);

	/* Build a light mesh to know where are the holes */
	vcn_mesh_t* mesh = alloca(vcn_mesh_get_memsize());
	vcn_mesh_init(mesh);
	vcn_mesh_get_simplest_from_model(mesh, model);
	
	/* Get holes and destroy mesh */
	uint32_t N_centroids;
	double* centroids =
		vcn_mesh_get_centroids_of_enveloped_areas(mesh, &N_centroids);
	vcn_mesh_finish(mesh);

	uint32_t N_holes = 0;
	double *holes = NULL;

	if (0 < N_centroids) {
		char* mask_centroids = NB_SOFT_MALLOC(N_centroids);
		memset(mask_centroids, 0, N_centroids);
		/* get mask */
		for (uint32_t i = 0; i < N_centroids; i++) {
			bool not_inside =
				!nb_msh3trg_is_vtx_inside(msh3trg,
							  &(centroids[i*2]));
			if (not_inside) {
				mask_centroids[i] = 1;
				N_holes += 1;
			}
		}

		if (0 < N_holes) {
			uint32_t i = 0;
			holes = malloc(2 * N_holes * sizeof(*holes));
			for (uint32_t j = 0; j < N_centroids; j++) {
				if (1 == mask_centroids[j]) {
					memcpy(&(holes[i*2]),
					       &(centroids[j*2]),
					       2 * sizeof(*holes));
					i += 1;
				}
			}
		}

		NB_SOFT_FREE(N_centroids, mask_centroids);
		free(centroids);
	}
	
	/* Build model with holes */
	model->H = N_holes;
	model->holes = holes;
}

vcn_model_t* vcn_model_create_from_msh3trg_with_disabled_trg
		(const vcn_msh3trg_t *const restrict msh3trg,
		 const bool *const restrict trg_enabled,
		 uint32_t *N_real_vtx_boundaries,
		 uint32_t **real_vtx_boundaries)
{
	/* Allocate model */
	vcn_model_t* model = vcn_model_create();
	/* Count segments and mark vertices used */
	model->M = 0;
	
	uint32_t vtx_memsize = msh3trg->N_vertices;
	char* vertices_used = NB_SOFT_MALLOC(vtx_memsize);
	memset(vertices_used, 0, vtx_memsize);

	char* vertices_bndr = NB_SOFT_MALLOC(vtx_memsize);
	memset(vertices_bndr, 0, vtx_memsize);

	for (uint32_t i = 0; i < msh3trg->N_triangles; i++){
		if (!trg_enabled[i])
			continue;
		uint32_t v1 = msh3trg->vertices_forming_triangles[i * 3];
		uint32_t v2 = msh3trg->vertices_forming_triangles[i*3+1];
		uint32_t v3 = msh3trg->vertices_forming_triangles[i*3+2];
		uint32_t nid = msh3trg->triangles_sharing_sides[i * 3];
		if (nid == msh3trg->N_triangles){
			model->M += 1;
			vertices_bndr[v1] = 1;
			vertices_bndr[v2] = 1;
			vertices_used[v1] = 1;
			vertices_used[v2] = 1;
		} else if (!trg_enabled[nid]) {
			model->M += 1;
			vertices_used[v1] = 1;
			vertices_used[v2] = 1;
		}

		nid = msh3trg->triangles_sharing_sides[i*3+1];
		if (nid == msh3trg->N_triangles) {
			model->M += 1;
			vertices_bndr[v2] = 1;
			vertices_bndr[v3] = 1;
			vertices_used[v2] = 1;
			vertices_used[v3] = 1;
		} else if (!trg_enabled[nid]) {
			model->M += 1;
			N_real_vtx_boundaries[0] += 1;
			vertices_used[v2] = 1;
			vertices_used[v3] = 1;
		}

		nid = msh3trg->triangles_sharing_sides[i*3+2];
		if (nid == msh3trg->N_triangles) {
			model->M += 1;
			vertices_bndr[v3] = 1;
			vertices_bndr[v1] = 1;
			vertices_used[v3] = 1;
			vertices_used[v1] = 1;
 		} else if (!trg_enabled[nid]) {
			model->M += 1;
			vertices_used[v3] = 1;
			vertices_used[v1] = 1;
		}
	}
	/* Count vertices */
	model->N = 0;
	for(uint32_t i = 0; i < msh3trg->N_vertices; i++){
		if(vertices_used[i]) 
			model->N += 1;
	}
	NB_SOFT_FREE(vtx_memsize, vertices_used);

	N_real_vtx_boundaries[0] = 0;
	for (uint32_t i = 0; i < msh3trg->N_vertices; i++) {
		if (vertices_bndr[i]) 
			N_real_vtx_boundaries[0] += 1;
	}

	/* Allocate segments and vertices */
	nb_model_alloc_vertices(model);
	nb_model_alloc_edges(model);

	real_vtx_boundaries[0] = malloc(N_real_vtx_boundaries[0] * 
					sizeof(uint32_t));

	/* Set vertices and segments */ 
	uint32_t idx_memsize = msh3trg->N_vertices * sizeof(uint32_t);
	uint32_t* vertices_idx = NB_SOFT_MALLOC(idx_memsize);
	memset(vertices_idx, 0, idx_memsize);

	for(uint32_t i = 0; i < msh3trg->N_vertices; i++)
		vertices_idx[i] = msh3trg->N_vertices;
  
	uint32_t real_vtx_cnt = 0;
	uint32_t vtx_counter = 0;
	uint32_t sgm_counter = 0;
	for(uint32_t i = 0; i < msh3trg->N_triangles; i++){
		if(!trg_enabled[i]) continue;
		uint32_t v1 = msh3trg->vertices_forming_triangles[i * 3];
		uint32_t v2 = msh3trg->vertices_forming_triangles[i*3+1];
		uint32_t v3 = msh3trg->vertices_forming_triangles[i*3+2];

		/* Check side 1 */
		uint32_t nid = msh3trg->triangles_sharing_sides[i * 3];
		bool include_side = false;
		if(nid == msh3trg->N_triangles)
			include_side = true;
		else if(!trg_enabled[nid])
			include_side = true;
		if(include_side){
			if(vertices_idx[v1] == msh3trg->N_vertices){
				memcpy(&(model->vertex[vtx_counter*2]),
				       &(msh3trg->vertices[v1*2]),
				       2* sizeof(*(model->vertex)));
				vertices_idx[v1] = vtx_counter++;
				if(vertices_bndr[v1])
					real_vtx_boundaries[0][real_vtx_cnt++] = 
						vertices_idx[v1];
			}
			if(vertices_idx[v2] == msh3trg->N_vertices){
				memcpy(&(model->vertex[vtx_counter*2]),
				       &(msh3trg->vertices[v2*2]),
				       2* sizeof(*(model->vertex)));
				vertices_idx[v2] = vtx_counter++;
				if(vertices_bndr[v2])
					real_vtx_boundaries[0][real_vtx_cnt++] =
						vertices_idx[v2];
			}
      
			model->edge[sgm_counter * 2] = vertices_idx[v1];
			model->edge[sgm_counter*2+1] = vertices_idx[v2];
			sgm_counter += 1;
		}

		/* Check side 2 */
		nid = msh3trg->triangles_sharing_sides[i*3+1];
		include_side = false;
		if(nid == msh3trg->N_triangles)
			include_side = true;
		else if(!trg_enabled[nid])
			include_side = true;
		if(include_side){
			if(vertices_idx[v2] == msh3trg->N_vertices){
				memcpy(&(model->vertex[vtx_counter*2]),
				       &(msh3trg->vertices[v2*2]),
				       2* sizeof(*(model->vertex)));
				vertices_idx[v2] = vtx_counter++;
				if(vertices_bndr[v2])
					real_vtx_boundaries[0][real_vtx_cnt++] =
						vertices_idx[v2];
			}
			if(vertices_idx[v3] == msh3trg->N_vertices){
				memcpy(&(model->vertex[vtx_counter*2]),
				       &(msh3trg->vertices[v3*2]),
				       2* sizeof(*(model->vertex)));
				vertices_idx[v3] = vtx_counter++;
				if(vertices_bndr[v3])
					real_vtx_boundaries[0][real_vtx_cnt++] = 
						vertices_idx[v3];
			}
      
			model->edge[sgm_counter * 2] = vertices_idx[v2];
			model->edge[sgm_counter*2+1] = vertices_idx[v3];
			sgm_counter += 1;
		}

		/* Check side 3 */
		nid = msh3trg->triangles_sharing_sides[i*3+2];
		include_side = false;
		if(nid == msh3trg->N_triangles)
			include_side = true;
		else if(!trg_enabled[nid])
			include_side = true;
		if(include_side){
			if(vertices_idx[v3] == msh3trg->N_vertices){
				memcpy(&(model->vertex[vtx_counter*2]),
				       &(msh3trg->vertices[v3*2]),
				       2* sizeof(*(model->vertex)));
				vertices_idx[v3] = vtx_counter++;
				if(vertices_bndr[v3])
					real_vtx_boundaries[0][real_vtx_cnt++] = 
						vertices_idx[v3];
			}
			if(vertices_idx[v1] == msh3trg->N_vertices){
				memcpy(&(model->vertex[vtx_counter*2]),
				       &(msh3trg->vertices[v1*2]),
				       2* sizeof(*(model->vertex)));
				vertices_idx[v1] = vtx_counter++;
				if(vertices_bndr[v1])
					real_vtx_boundaries[0][real_vtx_cnt++] =
						vertices_idx[v1];
			}
      
			model->edge[sgm_counter * 2] = vertices_idx[v3];
			model->edge[sgm_counter*2+1] = vertices_idx[v1];
			sgm_counter += 1;
		}
	}
	NB_SOFT_FREE(idx_memsize, vertices_idx);
	NB_SOFT_FREE(vtx_memsize, vertices_bndr);

	/* Return model */
	return model;
}

uint8_t vcn_model_save(const vcn_model_t *const model, const char* filename)
{
	/* Open file and verify if exist */
	FILE* fp = fopen(filename, "w");
	if (fp == NULL)
		return 1;
	/* Write number of vertices */
	fprintf(fp, "%u\n", model->N);
	/* Write vertices coordinates */
	for (uint32_t i = 0; i < model->N; i++)
		fprintf(fp, "%lf\t%lf\n", model->vertex[i * 2],
			model->vertex[i*2+1]);
  
	/* Write number of edge */
	fprintf(fp, "\n%u\n", model->M);
	/* Write segments */
	for (uint32_t i = 0; i<model->M; i++)
		fprintf(fp, "%u %u\n", model->edge[i * 2],
			model->edge[i*2+1]);

	/* Write number of holes */
	fprintf(fp, "\n%u\n", model->H);
	/* Write holes */
	for (uint32_t i = 0; i<model->H; i++)
		fprintf(fp, "%lf %lf\n", model->holes[i * 2],
			model->holes[i*2+1]);
	/* Close file */
	fclose(fp);
	/* Successful operation */
	return 0;
}

vcn_graph_t* vcn_model_get_vtx_graph(const vcn_model_t *const restrict model)
{
	vcn_graph_t *graph = vcn_graph_create();
  
	graph->N_adj = calloc(model->N, sizeof(*(graph->N_adj)));
	graph->adj = malloc(model->N * sizeof(*(graph->adj)));

	/* Compute connectivity matrix */
	for (uint32_t i = 0; i < model->M; i++) {
		uint32_t id_vi = model->edge[i * 2];
		uint32_t id_vj = model->edge[i*2+1];
		/* Count the number of connections */
		graph->N_adj[id_vi] += 1;
		graph->N_adj[id_vj] += 1;
	}
	for (uint32_t i = 0; i < model->N; i++)
		graph->adj[i] = malloc(graph->N_adj[i] *
				       sizeof(*(graph->adj[i])));

	/* Fill connectivity matrix */
	uint32_t rowc_memsize = model->N * sizeof(uint32_t);
	uint32_t* row_counter = NB_SOFT_MALLOC(rowc_memsize);
	memset(row_counter, 0, rowc_memsize);

	for (uint32_t i = 0; i < model->M; i++) {
		uint32_t id_vi = GET_1_EDGE_VTX(model, i);
		uint32_t id_vj = GET_2_EDGE_VTX(model, i);
		/* Fill matrix */
		graph->adj[id_vi][row_counter[id_vi]++] = id_vj;
		graph->adj[id_vj][row_counter[id_vj]++] = id_vi;
	}
	NB_SOFT_FREE(rowc_memsize, row_counter);

	return graph;
}

void vcn_model_set_enveloped_areas_as_holes(vcn_model_t* model)
{
	vcn_mesh_t* mesh = alloca(vcn_mesh_get_memsize());
	vcn_mesh_init(mesh);
	vcn_mesh_get_simplest_from_model(mesh, model);

	uint32_t N_holes;
	double* holes =
		vcn_mesh_get_centroids_of_enveloped_areas(mesh, &N_holes);
	vcn_mesh_finish(mesh);
	
	if (0 < model->H)
		free(model->holes);
	model->H = N_holes;
	model->holes = holes;
}

bool vcn_model_is_vtx_inside(const vcn_model_t *const model,
			     const double *const vtx)
{
	vcn_mesh_t* mesh = alloca(vcn_mesh_get_memsize());
	vcn_mesh_init(mesh);
	vcn_mesh_get_simplest_from_model(mesh, model);

	bool is_inside = vcn_mesh_is_vtx_inside(mesh, vtx);
	vcn_mesh_finish(mesh);
	return is_inside;
}

void vcn_model_get_enveloping_box(const vcn_model_t *const model,
				  double box[4])
{
	vcn_utils2D_get_enveloping_box(model->N, model->vertex,
				       2 * sizeof(*(model->vertex)),
				       vcn_utils2D_get_x_from_darray,
				       vcn_utils2D_get_y_from_darray,
				       box);
}

double* vcn_model_get_holes(const vcn_model_t *const model, 
			    uint32_t* N_holes /* Output */)
{
	double *holes;
	if (0 == model->H) {
		N_holes[0] = 0;
		holes = NULL;
	} else {
		N_holes[0] = model->H;
		holes = malloc(model->H * 2 * sizeof(*holes));
		memcpy(holes, model->holes, model->H * 2 * sizeof(*holes));
	}	
	return holes;
}

double* vcn_model_get_vertices(const vcn_model_t *const model, 
			       uint32_t* N_vertices /* Output */)
{
	double *vertices;
	if(0 == model->N) {
		N_vertices[0] = 0;
		vertices = NULL;
	} else {
		N_vertices[0] = model->N;
		vertices = malloc(2 * model->N * sizeof(*vertices));
		memcpy(vertices, model->vertex, 2 * model->N * sizeof(*vertices));
	}
	return vertices;
}

uint32_t vcn_model_get_vertex_id(const vcn_model_t *const model, double* vtx)
{
	uint32_t id = model->N;
	for (uint32_t i = 0; i < model->N; i++) {
		if (vcn_utils2D_get_dist2(vtx, &(model->vertex[i*2])) < NB_GEOMETRIC_TOL_POW2) {
			id = i;
			break;
		}
	}
	return id;
}

double* vcn_model_get_vertex_coordinate(const vcn_model_t *const model, uint32_t id)
{
	double *vtx = NULL;
	if (id < model->N) {
		vtx = malloc(2 * sizeof(*vtx));
		vtx[0] = model->vertex[id * 2];
		vtx[1] = model->vertex[id*2+1];
	}
	return vtx;
}

uint32_t vcn_model_get_edge_id(const vcn_model_t *const model, double* edge_vertices)
{
	uint32_t id = model->M;
	for (uint32_t i = 0; i < model->M; i++) {
		uint32_t id1 = model->edge[i * 2];
		uint32_t id2 = model->edge[i*2+1];
		if (vcn_utils2D_get_dist2(edge_vertices, 
				 &(model->vertex[id1 * 2])) <
		    NB_GEOMETRIC_TOL_POW2) {
			if (vcn_utils2D_get_dist2(&(edge_vertices[2]), 
					 &(model->vertex[id2*2])) <
			    NB_GEOMETRIC_TOL_POW2) {
				id = i;
				break;
			}
		} else if (vcn_utils2D_get_dist2(edge_vertices, 
					&(model->vertex[id2 * 2])) <
			   NB_GEOMETRIC_TOL_POW2){
			if (vcn_utils2D_get_dist2(&(edge_vertices[2]), 
					 &(model->vertex[id1*2])) <
			    NB_GEOMETRIC_TOL_POW2) {
				id = i;
				break;
			}
		}
	}
	return id;
}
double* vcn_model_get_edge_coordinates(const vcn_model_t *const model, uint32_t id)
{
	double *sgm = NULL;
	if(id < model->M) {
		uint32_t id1 = model->edge[id * 2];
		uint32_t id2 = model->edge[id*2+1];

		sgm = malloc(4 * sizeof(*sgm));
		sgm[0] = model->vertex[id1 * 2];
		sgm[1] = model->vertex[id1*2+1];
		sgm[2] = model->vertex[id2 * 2];
		sgm[3] = model->vertex[id2*2+1];
	}
	return sgm;
}

inline uint32_t vcn_model_get_number_of_vertices(const vcn_model_t *const model)
{
	return model->N;
}

inline uint32_t vcn_model_get_N_edges(const vcn_model_t *const model)
{
	return model->M;
}

inline uint32_t vcn_model_get_number_of_hole_seeds(const vcn_model_t *const model)
{
	return model->H;
}

double vcn_model_get_length_of_ith_edge(const vcn_model_t* model, uint32_t i)
{
	uint32_t v1 = model->edge[i * 2];
	uint32_t v2 = model->edge[i*2+1];
	return vcn_utils2D_get_dist(&(model->vertex[v1*2]),
			&(model->vertex[v2*2]));
}

double vcn_model_get_area(const vcn_model_t *const model)
{
	vcn_mesh_t* mesh = alloca(vcn_mesh_get_memsize());
	vcn_mesh_init(mesh);
	vcn_mesh_get_simplest_from_model(mesh, model);
	double area = vcn_mesh_get_area(mesh);
	vcn_mesh_finish(mesh);
	return area;
}

double vcn_model_get_sum_of_sgm_length(const vcn_model_t *const model)
{
	double sum = 0.0;
	for (uint32_t i = 0; i < model->M; i++) {
		uint32_t id1 = model->edge[i * 2];
		uint32_t id2 = model->edge[i*2+1];
		sum += vcn_utils2D_get_dist(&(model->vertex[id1 * 2]),
				&(model->vertex[id2 * 2]));
	}
	return sum;
}
