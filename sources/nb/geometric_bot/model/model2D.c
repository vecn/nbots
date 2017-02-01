#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>
#include <math.h>

#include "nb/math_bot.h"
#include "nb/memory_bot.h"
#include "nb/container_bot.h"
#include "nb/geometric_bot/utils2D.h"
#include "nb/geometric_bot/mesh/tessellator2D.h"
#include "nb/geometric_bot/mesh/modules2D/area_analizer.h"

#include "nb/geometric_bot/model/model2D_struct.h"
#include "nb/geometric_bot/model/model2D.h"

#define MAX(a, b) (((a)>(b))?(a):(b))

static inline void *nb_allocate_mem_model(void);
static void build_dynamic_graph(const nb_model_t *model,
				nb_membank_t *membank,
				nb_container_t **cnt_graph);
static void set_vtx_graph(nb_container_t **cnt_graph, nb_membank_t *membank,
			  nb_graph_t *graph);
static bool is_vtx_inside_polygon(const nb_model_t *const model,
				  const double vtx[2]);
static bool is_vtx_inside_mesh(const nb_model_t *const model,
			       const double vtx[2]);

uint16_t nb_model_get_memsize(void)
{
	return sizeof(nb_model_t);
}

void nb_model_init(void *model_ptr)
{
	memset(model_ptr, 0, nb_model_get_memsize());
}

void nb_model_copy(void *model_ptr, const void *src_model_ptr)
{
	nb_model_t *model = model_ptr;
	const nb_model_t *src_model = src_model_ptr;

	model->N = src_model->N;
	if (0 < model->N) {
		uint32_t size = 2 * model->N * sizeof(*(model->vertex));
		model->vertex = nb_allocate_mem(size);
		memcpy(model->vertex, src_model->vertex, size);
	}

	model->M = src_model->M;
	if (0 < model->M) {
		uint32_t size = 2 * model->M * sizeof(*(model->edge));
		model->edge = nb_allocate_mem(size);
		memcpy(model->edge, src_model->edge, size);
	}

	model->H = src_model->H;
	if (0 < model->H) {
		uint32_t size = 2 * model->H * sizeof(*(model->holes));
		model->holes = nb_allocate_mem(size);
		memcpy(model->holes, src_model->holes, size);
	}
}

void nb_model_finish(void *model_ptr)
{
	nb_model_clear(model_ptr);
}

void* nb_model_create(void)
{
	nb_model_t *model = nb_allocate_mem_model();
	nb_model_init(model);
	return model;
}

static inline void *nb_allocate_mem_model(void)
{
	uint32_t size = nb_model_get_memsize();
	return nb_allocate_mem(size);
}

void* nb_model_clone(const void *model_ptr)
{
	nb_model_t *model = nb_allocate_mem_model();
	nb_model_copy(model, model_ptr);
	return model;
}

void nb_model_destroy(void *model_ptr)
{
	nb_model_finish(model_ptr);
	nb_free_mem(model_ptr);
}

void nb_model_clear(void *model_ptr)
{
	nb_model_t *model = model_ptr;
	if (0 < model->N)
		nb_free_mem(model->vertex);
	model->N = 0;
	model->vertex = NULL;
	
	if (0 < model->M)
		nb_free_mem(model->edge);
	model->M = 0;
	model->edge = NULL;
	
	if (0 < model->H)
		nb_free_mem(model->holes);
	model->H = 0;
	model->holes = NULL;
}

nb_model_t* nb_model_load(const char* filename)
{
	/* OPPORTUNITY: Use custom format */
	nb_model_t* model = nb_model_create();
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

nb_model_t* nb_model_create_rectangle(double x_min, double y_min,
					double x_max, double y_max)
{
	nb_model_t* model = nb_model_create();
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

nb_model_t* nb_model_create_polygon(double radius,
				    double x_center,
				    double y_center,
				    uint32_t N_sides)
{
	nb_model_t *model = nb_model_create();
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

nb_model_t* nb_model_create_circle(double radius,
				   double x_center,
				   double y_center,
				   double side_length)
{
	double perimeter = 2.0 * NB_MATH_PI * radius;
	uint32_t n = (uint32_t) (perimeter / side_length + 0.5);
	if (10 > n)
		n = 10;
	return nb_model_create_polygon(radius, x_center, y_center, n);
}

uint8_t nb_model_save(const nb_model_t *const model, const char* filename)
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

void nb_model_load_vtx_graph(const nb_model_t *const model,
			      nb_graph_t *graph)
{
	uint32_t N = model->N;
	nb_container_type cnt_type = NB_SORTED;
	uint32_t cnt_size = nb_container_get_memsize(cnt_type);
	uint32_t memsize = N * (sizeof(void*) + cnt_size) 
		+ nb_membank_get_memsize();
	char *memblock = nb_soft_allocate_mem(memsize);

	nb_container_t **cnt_graph = (void*) memblock;
	for (uint32_t i = 0; i < N; i++) {
		cnt_graph[i] = (void*) (memblock + N * sizeof(void*) +
					i * cnt_size);
		nb_container_init(cnt_graph, cnt_type);
	}
	nb_membank_t *membank = (void*) (memblock + N * (sizeof(void*) +
							 cnt_size));
	nb_membank_init(membank, sizeof(uint32_t));

	build_dynamic_graph(model, membank, cnt_graph);

	graph->N = N;
	set_vtx_graph(cnt_graph, membank, graph);

	for (uint32_t i = 0; i < N; i++)
		nb_container_finish(cnt_graph);
	nb_membank_finish(membank);
	nb_soft_free_mem(memsize, memblock);
}

static void build_dynamic_graph(const nb_model_t *model,
				nb_membank_t *membank,
				nb_container_t **cnt_graph)
{
	for (uint32_t i = 0; i < model->M; i++) {
		uint32_t id1 = model->edge[i * 2];
		uint32_t id2 = model->edge[i*2+1];

		uint32_t *aux1 = nb_membank_allocate_mem(membank);
		*aux1 = id1;
		uint32_t *aux2 = nb_membank_allocate_mem(membank);
		*aux2 = id2;

		nb_container_insert(cnt_graph[id1], aux2);
		nb_container_insert(cnt_graph[id2], aux1);
	}
}

static void set_vtx_graph(nb_container_t **cnt_graph, nb_membank_t *membank,
			       nb_graph_t *graph)
{
	uint32_t N_total_adj = 0;
	for (uint32_t i = 0; i < graph->N; i++)
		N_total_adj += nb_container_get_length(cnt_graph[i]);
	uint32_t memsize = graph->N * (sizeof(*(graph->N_adj)) +
				       sizeof(*(graph->adj))) +
		N_total_adj * sizeof(**(graph->adj));

	char *memblock = nb_allocate_mem(memsize);

	graph->N_adj = (void*) memblock;
	graph->adj = (void*) (memblock + sizeof(*(graph->N_adj)));

	memblock += graph->N * (sizeof(*(graph->N_adj)) +
				sizeof(*(graph->adj)));
	for (uint32_t i = 0; i < graph->N; i++) {
		uint32_t N_adj = nb_container_get_length(cnt_graph[i]);
		graph->N_adj[i] = N_adj;
		graph->adj[i] = (void*) memblock;
		memblock += N_adj * sizeof(**(graph->adj));

		uint32_t j = 0;
		while (nb_container_is_not_empty(cnt_graph[i])) {
			uint32_t *id = nb_container_delete_first(cnt_graph[i]);
			graph->adj[i][j] = *id;
			j += 1;
			nb_membank_free_mem(membank, id);
		}
	}
}

void nb_model_set_enveloped_areas_as_holes(nb_model_t* model)
{
	nb_tessellator2D_t* mesh = nb_allocate_on_stack(nb_tessellator2D_get_memsize());
	nb_tessellator2D_init(mesh);
	nb_tessellator2D_get_simplest_from_model(mesh, model);

	uint32_t N_holes;
	double* holes =
		nb_tessellator2D_get_centroids_of_enveloped_areas(mesh, &N_holes);
	nb_tessellator2D_finish(mesh);
	
	if (0 < model->H)
		nb_free_mem(model->holes);

	model->H = N_holes;
	model->holes = holes;
}

bool nb_model_is_vtx_inside(const nb_model_t *const model,
			     const double *const vtx)
{
	bool out;
	if (nb_model_get_N_edges(model) < 1)/* AQUI VOY */
		out = is_vtx_inside_polygon(model, vtx);
	else
		out = is_vtx_inside_mesh(model, vtx);
	return out;
}

static bool is_vtx_inside_polygon(const nb_model_t *const model,
				  const double vtx[2])
{
	uint32_t N = nb_model_get_number_of_vertices(model);
	return nb_utils2D_pnt_lies_in_poly(N, model->vertex, vtx);
}

static bool is_vtx_inside_mesh(const nb_model_t *const model,
			       const double vtx[2])
{
	nb_tessellator2D_t* mesh = nb_allocate_on_stack(nb_tessellator2D_get_memsize());
	nb_tessellator2D_init(mesh);
	nb_tessellator2D_get_simplest_from_model(mesh, model);

	bool is_inside = nb_tessellator2D_is_vtx_inside(mesh, vtx);
	nb_tessellator2D_finish(mesh);
	return is_inside;
}

void nb_model_get_enveloping_box(const nb_model_t *const model,
				  double box[4])
{
	nb_utils2D_get_enveloping_box(model->N, model->vertex,
				       2 * sizeof(*(model->vertex)),
				       nb_utils2D_get_x_from_darray,
				       nb_utils2D_get_y_from_darray,
				       box);
}

double* nb_model_get_holes(const nb_model_t *const model, 
			    uint32_t* N_holes /* Output */)
{
	double *holes;
	if (0 == model->H) {
		N_holes[0] = 0;
		holes = NULL;
	} else {
		N_holes[0] = model->H;
		holes = nb_allocate_mem(model->H * 2 * sizeof(*holes));
		memcpy(holes, model->holes, model->H * 2 * sizeof(*holes));
	}	
	return holes;
}

double* nb_model_get_vertices(const nb_model_t *const model, 
			       uint32_t* N_vertices /* Output */)
{
	double *vertices;
	if(0 == model->N) {
		N_vertices[0] = 0;
		vertices = NULL;
	} else {
		N_vertices[0] = model->N;
		vertices = nb_allocate_mem(2 * model->N * sizeof(*vertices));
		memcpy(vertices, model->vertex, 2 * model->N * sizeof(*vertices));
	}
	return vertices;
}

uint32_t nb_model_get_vertex_id(const nb_model_t *const model, double* vtx)
{
	uint32_t id = model->N;
	for (uint32_t i = 0; i < model->N; i++) {
		if (nb_utils2D_get_dist2(vtx, &(model->vertex[i*2])) < NB_GEOMETRIC_TOL_POW2) {
		  id = i;
			break;
		}
	}
	return id;
}

void nb_model_get_vertex_coordinate(const nb_model_t *const model, uint32_t id,
				    double vtx[2])
{
	if (id < model->N) {
		vtx[0] = model->vertex[id * 2];
		vtx[1] = model->vertex[id*2+1];
	}
}

uint32_t nb_model_get_edge_id(const nb_model_t *const model, double* edge_vertices)
{
	uint32_t id = model->M;
	for (uint32_t i = 0; i < model->M; i++) {
		uint32_t id1 = model->edge[i * 2];
		uint32_t id2 = model->edge[i*2+1];
		if (nb_utils2D_get_dist2(edge_vertices, 
				 &(model->vertex[id1 * 2])) <
		    NB_GEOMETRIC_TOL_POW2) {
			if (nb_utils2D_get_dist2(&(edge_vertices[2]), 
					 &(model->vertex[id2*2])) <
			    NB_GEOMETRIC_TOL_POW2) {
				id = i;
				break;
			}
		} else if (nb_utils2D_get_dist2(edge_vertices, 
					&(model->vertex[id2 * 2])) <
			   NB_GEOMETRIC_TOL_POW2){
			if (nb_utils2D_get_dist2(&(edge_vertices[2]), 
					 &(model->vertex[id1*2])) <
			    NB_GEOMETRIC_TOL_POW2) {
				id = i;
				break;
			}
		}
	}
	return id;
}
void nb_model_get_edge_coordinates(const nb_model_t *const model, uint32_t id,
				   double v1[2], double v2[2])
{
	if(id < model->M) {
		uint32_t id1 = model->edge[id * 2];
		uint32_t id2 = model->edge[id*2+1];

		v1[0] = model->vertex[id1 * 2];
		v1[1] = model->vertex[id1*2+1];
		v2[0] = model->vertex[id2 * 2];
		v2[1] = model->vertex[id2*2+1];
	}
}

uint32_t nb_model_get_number_of_vertices(const nb_model_t *const model)
{
	return model->N;
}

uint32_t nb_model_get_N_edges(const nb_model_t *const model)
{
	return model->M;
}

uint32_t nb_model_get_number_of_hole_seeds(const nb_model_t *const model)
{
	return model->H;
}

double nb_model_get_length_of_ith_edge(const nb_model_t* model, uint32_t i)
{
	uint32_t v1 = model->edge[i * 2];
	uint32_t v2 = model->edge[i*2+1];
	return nb_utils2D_get_dist(&(model->vertex[v1*2]),
			&(model->vertex[v2*2]));
}

double nb_model_get_area(const nb_model_t *const model)
{
	nb_tessellator2D_t* mesh = nb_allocate_on_stack(nb_tessellator2D_get_memsize());
	nb_tessellator2D_init(mesh);
	nb_tessellator2D_get_simplest_from_model(mesh, model);
	double area = nb_tessellator2D_get_area(mesh);
	nb_tessellator2D_finish(mesh);
	return area;
}

double nb_model_get_sum_of_sgm_length(const nb_model_t *const model)
{
	double sum = 0.0;
	for (uint32_t i = 0; i < model->M; i++) {
		uint32_t id1 = model->edge[i * 2];
		uint32_t id2 = model->edge[i*2+1];
		sum += nb_utils2D_get_dist(&(model->vertex[id1 * 2]),
				&(model->vertex[id2 * 2]));
	}
	return sum;
}

double* nb_model_get_centroids_of_enveloped_areas(const nb_model_t *model,
						  uint32_t *N_centroids)
{
	uint32_t memsize = nb_tessellator2D_get_memsize();
	nb_tessellator2D_t* mesh = nb_allocate_on_stack(memsize);
	nb_tessellator2D_init(mesh);
	nb_tessellator2D_get_simplest_from_model(mesh, model);
	
	double* centroids =
		nb_tessellator2D_get_centroids_of_enveloped_areas(mesh, N_centroids);
	nb_tessellator2D_finish(mesh);
	return centroids;
}
