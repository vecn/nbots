#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include "nb/memory_bot.h"
#include "nb/container_bot.h"
#include "nb/geometric_bot/point2D.h"
#include "nb/geometric_bot/utils2D.h"
#include "nb/geometric_bot/knn/bins2D.h"
#include "nb/geometric_bot/mesh/mesh2D.h"
#include "nb/geometric_bot/mesh/dewall.h"

#include "mesh2D_structs.h"

#define POW2(a) ((a)*(a))
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

typedef struct {
	double min;
	double max;
} interval_t;

typedef struct {
	bool using_bins;
	uint32_t N;
	nb_bins2D_t *bins;
	msh_vtx_t** vtx_array;
	int8_t axe;
	double alpha;
	int N_half;
} search_vtx_t;

static  void mesh_init_disp_and_scale
                            (nb_tessellator2D__t *const mesh,
			     uint32_t N_vertices,
			     const double *const vertices);
static int8_t compare_using_axe(const void *const v1_ptr,
				const void *const v2_ptr,
				const void *const axe_ptr);
static msh_vtx_t* get_1st_vtx(search_vtx_t *search_vtx);
static msh_vtx_t* get_2nd_vtx(search_vtx_t* search_vtx, msh_vtx_t* v1);
static msh_vtx_t* get_2nd_vtx_exahustive_search(uint32_t N, 
						msh_vtx_t** vertices,
						msh_vtx_t* v1);
static bool compare_2nd_vtx(const nb_point2D_t *const p_ref,
			    const nb_point2D_t *const p,
			    const void *const data);
static msh_vtx_t* get_2nd_vtx_using_bins(nb_bins2D_t* bins, msh_vtx_t* v1,
					 int8_t axe, double alpha);
static double swap_vtx_if_minimize_delaunay_distance
                                      (const nb_container_t *const edges,
				       const msh_vtx_t *const v1,
				       const msh_vtx_t *const v2,
				       msh_vtx_t** v3, double min_dist,
				       msh_vtx_t* vtx, double dist);
static double handle_cocircularity(const nb_container_t *const edges,
				   const msh_vtx_t *const v1,
				   const msh_vtx_t *const v2,
				   msh_vtx_t** v3, double min_dist,
				   msh_vtx_t *vtx, double dist);
static bool proposed_trg_intersects_edge(const msh_vtx_t *const v1,
					 const msh_vtx_t *const v2,
					 const msh_vtx_t *const v3,
					 const nb_container_t *const edges);
static msh_vtx_t* get_3rd_vtx(const nb_container_t *const edges,
			      const search_vtx_t *const search_vtx,
			      const msh_vtx_t *const  v1,
			      const msh_vtx_t *const  v2);
static msh_vtx_t* get_3rd_vtx_exahustive_search(const nb_container_t *const edges,
						uint32_t N, 
						msh_vtx_t** vertices,
						const msh_vtx_t *const v1,
						const msh_vtx_t *const v2);
static double set_v3(const nb_container_t *const edges,
		     const msh_vtx_t *const v1, const msh_vtx_t *const v2,
		     msh_vtx_t* v3_candidate, msh_vtx_t** v3, double min_dist);
static msh_vtx_t* get_3rd_vtx_using_bins(const nb_container_t *const edges,
					 const nb_bins2D_t *const bins,
					 const msh_vtx_t *const  v1,
					 const msh_vtx_t *const  v2);
static msh_trg_t* create_1st_trg(nb_tessellator2D__t *mesh, search_vtx_t* search_vtx);
static msh_trg_t* create_trg(nb_tessellator2D__t *mesh,
			     const search_vtx_t *search_vtx,
			     msh_edge_t *edge);
static void update_AFL(nb_container_t *AFL,
		       const msh_edge_t *const edge);
static uint32_t dewall_recursion
                       (nb_tessellator2D__t* mesh, uint32_t N,
			msh_vtx_t** vertices,
			uint16_t deep_level,
			nb_container_t* AFL);
static uint32_t triangulate_wall(nb_tessellator2D__t *mesh,
				 const search_vtx_t *search_vtx,
				 nb_container_t *AFL_alpha,
				 nb_container_t *AFL_1,
				 nb_container_t *AFL_2);
static void update_AFLs(const msh_trg_t *const trg,
			const msh_edge_t *const edge,
			nb_container_t *AFL_alpha,
			nb_container_t *AFL_1,
			nb_container_t *AFL_2,
			const search_vtx_t *search_vtx);
static double get_alpha(msh_vtx_t **vertices,
			int8_t axe, uint32_t N_mid);
static nb_container_t* get_side_AFL(const msh_edge_t *const edge,
				    const search_vtx_t *search_vtx,
				    nb_container_t *AFL_neg,
				    nb_container_t *AFL_zero,
				    nb_container_t *AFL_pos);

static uint32_t split_vtx_array(const nb_container_t *const AFL,
				uint32_t N, msh_vtx_t **vertices,
				int8_t axe);
static void set_interval(interval_t *interval,
			 const nb_container_t *const AFL,
			 uint32_t N, msh_vtx_t **vertices,
			 int8_t axe);
static uint32_t get_mid_inside_interval(interval_t *interval, uint32_t N,
					msh_vtx_t **vertices, int8_t axe);
static bool not_space_for_alpha(msh_vtx_t **vertices, uint32_t i, int8_t axe);
static void init_search_vtx(search_vtx_t *search_vtx, uint32_t N,
			    msh_vtx_t **vertices, int8_t axe,
			    double alpha, uint32_t N_half);
static bool set_first_trg_into_AFL(nb_tessellator2D__t *mesh,
				   search_vtx_t *search_vtx,
				   nb_container_t *AFL);
static void clear_search_vtx(search_vtx_t *search_vtx);
static uint32_t dewall(nb_tessellator2D__t* mesh);

void nb_tessellator2D__get_delaunay(nb_tessellator2D__t *mesh, uint32_t N_vertices,
			  const double *const vertices)
{
	nb_tessellator2D__clear(mesh);
	if (0 == N_vertices)
		goto EXIT;

	mesh_init_disp_and_scale(mesh, N_vertices, vertices);
	mesh->N_input_vtx = N_vertices;
	mesh->input_vtx =
		nb_allocate_zero_mem(mesh->N_input_vtx *
				     sizeof(*(mesh->input_vtx)));

	for (uint32_t i = 0; i < N_vertices; i++) {
		msh_vtx_t* vtx = mvtx_create(mesh);
		mesh->input_vtx[i] = vtx;
		mvtx_set_type_origin(vtx, INPUT);
		vtx->x[0] = mesh->scale *
			(vertices[i * 2] - mesh->xdisp);
		vtx->x[1] = mesh->scale *
			(vertices[i*2+1] - mesh->ydisp);
		nb_bins2D_insert(mesh->ug_vtx, vtx);
	}
	dewall(mesh);
EXIT:
	return;
}

static inline void mesh_init_disp_and_scale
                            (nb_tessellator2D__t *const restrict mesh,
			     uint32_t N_vertices,
			     const double *const restrict vertices)
{
	double box[4];
	nb_utils2D_get_enveloping_box(N_vertices, vertices,
				       2 * sizeof(*vertices),
				       nb_utils2D_get_x_from_darray,
				       nb_utils2D_get_y_from_darray,
				       box);

	mesh->xdisp = (box[0] + box[2]) / 2.0;
	mesh->ydisp = (box[1] + box[3]) / 2.0;
	mesh->scale = 100.0 / MAX(box[2] - box[0], box[3] - box[1]);
}

static inline int8_t compare_using_axe(const void *const v1_ptr,
				       const void *const v2_ptr,
				       const void *const axe_ptr)
{
	const msh_vtx_t *restrict v1 = *((msh_vtx_t**)v1_ptr);
	const msh_vtx_t *restrict v2 = *((msh_vtx_t**)v2_ptr);
	int8_t axe = *((int8_t*)axe_ptr);
	int8_t out;
	if (v1->x[axe] < v2->x[axe])
		out = -1;
	else if (v1->x[axe] > v2->x[axe])
		out = 1;
	else
		out = 0;
	return out;
}

static inline msh_vtx_t* get_2nd_vtx(search_vtx_t* search_vtx, msh_vtx_t* v1)
{
	msh_vtx_t *v2;
	if (search_vtx->using_bins)
		v2 = get_2nd_vtx_using_bins(search_vtx->bins, v1,
					    search_vtx->axe,
					    search_vtx->alpha);
	else
		v2 = get_2nd_vtx_exahustive_search(search_vtx->N_half,
						   search_vtx->vtx_array, v1);
	return v2;
}

static inline msh_vtx_t* get_2nd_vtx_exahustive_search
                                               (uint32_t N, 
						msh_vtx_t** vertices,
						msh_vtx_t* v1)
{
	double min_dist = nb_utils2D_get_dist2(v1->x, vertices[0]->x);
	uint32_t id =  0;
	for (uint32_t i = 1; i < N; i++) {
		double dist = nb_utils2D_get_dist2(v1->x, vertices[i]->x);
		if (dist < min_dist) {
			min_dist = dist;
			id = i;
		}
	}
	return vertices[id];
}

static inline bool compare_2nd_vtx(const nb_point2D_t *const restrict p_ref,
				   const nb_point2D_t *const restrict p,
				   const void *const restrict data)
{
	register double alpha = ((double*)data)[0];
	register int8_t axe = (uint8_t)((double*)data)[1];
	return (alpha - p->x[axe] > 0.0);

}

static inline msh_vtx_t* get_2nd_vtx_using_bins(nb_bins2D_t* bins,
						msh_vtx_t* v1, int8_t axe,
						double alpha)
{
	double filter_data[2];
	filter_data[0] = alpha;
	filter_data[1] = (double) axe;

	msh_vtx_t* v2 = NULL;
	double dist;
	nb_bins2D_set_filter(bins, compare_2nd_vtx);
	nb_bins2D_set_filter_data(bins, filter_data);
	nb_bins2D_get_knn(bins, v1, 1, &v2, &dist);
	return v2;
}

static double swap_vtx_if_minimize_delaunay_distance
                                      (const nb_container_t *const edges,
				       const msh_vtx_t *const v1,
				       const msh_vtx_t *const v2,
				       msh_vtx_t** v3, double min_dist,
				       msh_vtx_t *vtx, double dist)
{
	if (fabs(dist - min_dist) < NB_GEOMETRIC_TOL) {
		min_dist = handle_cocircularity(edges, v1, v2, v3,
						min_dist, vtx, dist);
	} else {
		if (dist < min_dist) {
			*v3 = vtx;
			min_dist = dist;    
		}
	}
	return min_dist;
}

static double handle_cocircularity(const nb_container_t *const edges,
				   const msh_vtx_t *const v1,
				   const msh_vtx_t *const v2,
				   msh_vtx_t** v3, double min_dist,
				   msh_vtx_t *vtx, double dist)
{
	bool v3_intersect = proposed_trg_intersects_edge(v1, v2, *v3, edges);
	if (v3_intersect) {
		*v3 = vtx;
		min_dist = dist;
	}
	return min_dist;
}

static bool proposed_trg_intersects_edge(const msh_vtx_t *const v1,
					 const msh_vtx_t *const v2,
					 const msh_vtx_t *const v3,
					 const nb_container_t *const edges)
{
	bool intersects = false;
	uint32_t iter_size = nb_iterator_get_memsize();
	nb_iterator_t *iter = nb_soft_allocate_mem(iter_size);
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, edges);
	while (nb_iterator_has_more(iter)) {
		const msh_edge_t *edge = nb_iterator_get_next(iter);
		bool intersected =
			NB_INTERSECTED == nb_utils2D_get_sgm_intersection
							(v1->x, v3->x,
							 edge->v1->x,
							 edge->v2->x,
							 NULL);
		if (intersected) {
			intersects = true;
			break;
		} else {
			intersected = (NB_INTERSECTED == 
				       nb_utils2D_get_sgm_intersection
				       			(v2->x, v3->x,
							 edge->v1->x,
							 edge->v2->x,
							 NULL));
			if (intersected) {
				intersects = true;
				break;
			}
		}
	}
	nb_iterator_finish(iter);
	nb_soft_free_mem(iter_size, iter);
	return intersects;
}

static inline msh_vtx_t* get_3rd_vtx(const nb_container_t *const edges,
				     const search_vtx_t *const search_vtx,
				     const msh_vtx_t *const  v1,
				     const msh_vtx_t *const  v2)
{
	msh_vtx_t *v3;
	if (search_vtx->using_bins)
		v3 = get_3rd_vtx_using_bins(edges, search_vtx->bins, v1, v2);
	else
		v3 = get_3rd_vtx_exahustive_search(edges, search_vtx->N,
						   search_vtx->vtx_array,
						   v1, v2);
	return v3;
}

static msh_vtx_t* get_3rd_vtx_exahustive_search(const nb_container_t *const edges,
						uint32_t N, 
						msh_vtx_t** vertices,
						const msh_vtx_t *const v1,
						const msh_vtx_t *const v2)
{
	msh_vtx_t *v3 = NULL;
	double min_dist = 0.0;
	for (uint32_t i = 0; i < N; i++) {
		if (vertices[i] != v1 && vertices[i] != v2) {
			if (nb_utils2D_is_in_half_side(v1->x, v2->x,
							vertices[i]->x))
				min_dist = set_v3(edges, v1, v2, vertices[i],
						  &v3, min_dist);
		}
	}
	return v3;
}

static double set_v3(const nb_container_t *const edges,
		     const msh_vtx_t *const v1, const msh_vtx_t *const v2,
		     msh_vtx_t* v3_candidate, msh_vtx_t** v3, double min_dist)
{
	double dist = nb_utils2D_get_delaunay_dist(v1->x, v2->x,
						    v3_candidate->x);
	if (NULL == *v3) {
		/* Initialize minimum distance */
		min_dist = dist;
		*v3 = v3_candidate;
	} else {
		min_dist = 
			swap_vtx_if_minimize_delaunay_distance(edges,
							       v1, v2, v3,
							       min_dist,
							       v3_candidate,
							       dist);
	}
	return min_dist;

}

static msh_vtx_t* get_3rd_vtx_using_bins(const nb_container_t *const edges,
					 const nb_bins2D_t *const restrict bins,
					 const msh_vtx_t *const restrict v1,
					 const msh_vtx_t *const restrict v2)
{
	uint32_t cnt_size = nb_container_get_memsize(NB_QUEUE);
	nb_container_t* vertices = nb_soft_allocate_mem(cnt_size);
	nb_container_init(vertices, NB_QUEUE);

	nb_bins2D_get_candidate_points_to_min_delaunay(bins, v1, v2,
							vertices);
  	msh_vtx_t *v3 = NULL;
	double min_dist = 0.0;
	while (nb_container_is_not_empty(vertices)) {
		msh_vtx_t* vtx = nb_container_delete_first(vertices);
		min_dist = set_v3(edges, v1, v2, vtx, &v3, min_dist);
	}
	nb_container_finish(vertices);
	nb_soft_free_mem(cnt_size, vertices);
       	return v3;
}

static msh_trg_t* create_1st_trg(nb_tessellator2D__t *mesh, search_vtx_t* search_vtx)
{
	msh_vtx_t *v1 = get_1st_vtx(search_vtx);
	msh_vtx_t *v2 = get_2nd_vtx(search_vtx, v1);
	msh_vtx_t *v3 = get_3rd_vtx(mesh->ht_edge, search_vtx, v1, v2);
  
	if (NULL == v3) {
		/* Check in the other halfspace */
		msh_vtx_t *aux_vtx = v1;
		v1 = v2;
		v2 = aux_vtx;
		v3 = get_3rd_vtx(mesh->ht_edge, search_vtx, v1, v2);
	}

	msh_trg_t *trg;
	if (NULL != v3) {
		trg = mtrg_allocate_zero_mem(mesh);
		trg->v1 = v1;
		trg->v2 = v2;
		trg->v3 = v3;
	} else {
		/* All points are collinear */
		trg = NULL;
	}
	return trg;
}

static inline msh_vtx_t* get_1st_vtx(search_vtx_t *search_vtx)
{
	/* Closest vertex to the wall */
	return search_vtx->vtx_array[search_vtx->N_half];
}

static msh_trg_t* create_trg(nb_tessellator2D__t *mesh,
			     const search_vtx_t *search_vtx,
			     msh_edge_t *edge)
{
	/* Select the correct segment orientation */
	msh_vtx_t *restrict v1;
	msh_vtx_t *restrict v2;
	if (NULL == edge->t1) {
		v1 = edge->v1;
		v2 = edge->v2;
	} else {
		v1 = edge->v2;
		v2 = edge->v1;
	}

	msh_vtx_t *restrict v3 = get_3rd_vtx(mesh->ht_edge,
					     search_vtx, v1, v2);

	msh_trg_t *trg = NULL;
	if(NULL != v3) {
		trg = mtrg_allocate_zero_mem(mesh);
		trg->v1 = v1;
		trg->v2 = v2;
		trg->v3 = v3;
	}
	return trg;
}

static inline void update_AFL(nb_container_t *AFL,
			      const msh_edge_t *const edge)
{
	if (NULL == nb_container_exist(AFL, edge))
		nb_container_insert(AFL, edge);
	else
		nb_container_delete(AFL, edge);
}

static uint32_t dewall_recursion
                       (nb_tessellator2D__t* mesh, uint32_t N,
			msh_vtx_t** vertices,
			uint16_t deep_level,
			nb_container_t* AFL)
{
	uint32_t N_trg = 0;
	if (N < 3)
		goto EXIT;

	/* Ascending sorting of vertices depending on the axe divisor */
	int8_t axe = deep_level % 2; /* 0:X, 1:Y */

	nb_qsort_wd(vertices, N, sizeof(*vertices),
		     compare_using_axe, &axe);
	uint32_t N_mid = split_vtx_array(AFL, N, vertices, axe);
	double alpha;
	if (N > 3) {
		alpha = get_alpha(vertices, axe, N_mid);
	} else {
		double min = vertices[0]->x[axe];
		double max = vertices[0]->x[axe];
		if (min > vertices[1]->x[axe])
			min = vertices[1]->x[axe];
		if (max < vertices[1]->x[axe])
			max = vertices[1]->x[axe];
		if (max < vertices[2]->x[axe])
			max = vertices[2]->x[axe];
		alpha = (min + max)/2.0;
	}

	search_vtx_t search_vtx;
	init_search_vtx(&search_vtx, N, vertices, axe, alpha, N_mid);

	/* Initialize triangles counter */
	uint32_t n_trg_alpha = 0;

	/* Create first triangle */
	if (nb_container_is_empty(AFL)) {
		if (set_first_trg_into_AFL(mesh, &search_vtx, AFL))
			n_trg_alpha += 1;
		else
			goto EXIT;
	}

	/* Initialize Action face lists */
	uint32_t AFL_size = nb_container_get_memsize(NB_HASH);
	uint32_t memsize = 3 * AFL_size;
	char *memblock = nb_allocate_mem(memsize);

	nb_container_t* AFL_alpha = (void*) memblock;
	nb_container_init(AFL_alpha, NB_HASH);
	nb_container_set_key_generator(AFL_alpha, hash_key_edge);

	nb_container_t* AFL_1 = (void*) (memblock + AFL_size);
	nb_container_init(AFL_1, NB_HASH);
	nb_container_set_key_generator(AFL_1, hash_key_edge);

	nb_container_t* AFL_2 = (void*) (memblock + 2 * AFL_size);
	nb_container_init(AFL_2, NB_HASH);
	nb_container_set_key_generator(AFL_2, hash_key_edge);

	/* Redistribute segments from the main AFL to the three AFLs */
	while (nb_container_is_not_empty(AFL)) {
		msh_edge_t* edge = nb_container_delete_first(AFL);
		nb_container_t *side_AFL = get_side_AFL(edge, &search_vtx, 
							AFL_1, AFL_alpha,
							AFL_2);
		nb_container_insert(side_AFL, edge);
	}
	
	n_trg_alpha += triangulate_wall(mesh, &search_vtx, AFL_alpha,
					AFL_1, AFL_2);
	nb_container_finish(AFL_alpha);

	clear_search_vtx(&search_vtx);

	/* Recursive triangulation */
	uint32_t n_trg_1 = 0;
	if (nb_container_is_not_empty(AFL_1))
		n_trg_1 += dewall_recursion(mesh, N_mid, vertices,
					    deep_level + 1, AFL_1);
	nb_container_finish(AFL_1);

	uint32_t n_trg_2 = 0;
	if (nb_container_is_not_empty(AFL_2))
		n_trg_2 += dewall_recursion(mesh, N - N_mid, &(vertices[N_mid]),
					    deep_level + 1, AFL_2);
	nb_container_finish(AFL_2);

	N_trg = n_trg_alpha + n_trg_1 + n_trg_2;

	nb_free_mem(memblock);
EXIT:
	return N_trg;
}

static inline double get_alpha(msh_vtx_t **vertices,
			       int8_t axe, uint32_t N_mid)
{
	double alpha;
	if (0 == N_mid)
		alpha = vertices[0]->x[axe];
	else
		alpha = 0.6 * vertices[N_mid]->x[axe] +
			0.4 * vertices[N_mid-1]->x[axe];
	return alpha;
}

static inline nb_container_t* get_side_AFL(const msh_edge_t *const edge,
					   const search_vtx_t *search_vtx,
					   nb_container_t *AFL_neg,
					   nb_container_t *AFL_zero,
					   nb_container_t *AFL_pos)
{
	bool v1 = (edge->v1->x[search_vtx->axe] - search_vtx->alpha > 0.0);
	bool v2 = (edge->v2->x[search_vtx->axe] - search_vtx->alpha > 0.0);
	nb_container_t *AFL;
	if (v1 != v2)
		AFL = AFL_zero;
	else if (v1)
		AFL = AFL_pos;
	else 
		AFL = AFL_neg;
	return AFL;
}

static uint32_t split_vtx_array(const nb_container_t *const AFL,
				uint32_t N, msh_vtx_t **vertices,
				int8_t axe)
{
	interval_t interval;
	set_interval(&interval, AFL, N, vertices, axe);

	uint32_t mid_in_interval =
		get_mid_inside_interval(&interval, N, vertices, axe);

	uint32_t N_mid = mid_in_interval;
	while (not_space_for_alpha(vertices, N_mid, axe)) {
		N_mid += 1;
		if (N == N_mid)
			break;
	}
	if (N == N_mid) {
		N_mid = mid_in_interval;
		while (not_space_for_alpha(vertices, N_mid, axe)) {
			N_mid -= 1;
			if (0 == N_mid)
				break;
		}
	}
	return N_mid;
}

static void set_interval(interval_t *interval,
			 const nb_container_t *const AFL,
			 uint32_t N, msh_vtx_t **vertices,
			 int8_t axe)
{
	if (nb_container_is_not_empty(AFL)) {
		uint32_t iter_size = nb_iterator_get_memsize();
		nb_iterator_t *iter = nb_soft_allocate_mem(iter_size);
		nb_iterator_init(iter);
		nb_iterator_set_container(iter, AFL);
		const msh_edge_t *edge = nb_iterator_get_next(iter);
		if (edge->v1->x[axe] < edge->v2->x[axe]) {
			interval->min = edge->v1->x[axe];
			interval->max = edge->v2->x[axe];
		} else {
			interval->min = edge->v2->x[axe];
			interval->max = edge->v1->x[axe];
		}
		while (nb_iterator_has_more(iter)) {
			edge = nb_iterator_get_next(iter);
			if (edge->v1->x[axe] < edge->v2->x[axe]) {
				interval->min = MIN(interval->min,
						    edge->v1->x[axe]);
				interval->max = MAX(interval->max,
						    edge->v2->x[axe]);
			} else {
				interval->min = MIN(interval->min,
						    edge->v2->x[axe]);
				interval->max = MAX(interval->max,
						    edge->v1->x[axe]);
			}
		}
		nb_iterator_finish(iter);
		nb_soft_free_mem(iter_size, iter);
	} else {
		interval->min = vertices[0]->x[axe];
		interval->max = vertices[N-1]->x[axe];
	}
}

static uint32_t get_mid_inside_interval(interval_t *interval, uint32_t N,
					msh_vtx_t **vertices, int8_t axe)
{
	uint32_t N_mid = N/2;
	while (vertices[N_mid]->x[axe] < interval->min) {
		N_mid += 1;
		if (N == N_mid)
			break;
	}
	while (vertices[N_mid]->x[axe] > interval->max) {
		N_mid -= 1;
		if (0 == N_mid)
			break;
	}
	return N_mid;
}

static inline bool not_space_for_alpha(msh_vtx_t **vertices, uint32_t i,
				       int8_t axe)
{
	return (vertices[i]->x[axe] - vertices[i-1]->x[axe])
		< NB_GEOMETRIC_TOL;
}

static void init_search_vtx(search_vtx_t *search_vtx, uint32_t N,
			    msh_vtx_t **vertices, int8_t axe,
			    double alpha, uint32_t N_half)
{
	search_vtx->N = N;
	search_vtx->vtx_array = vertices;
	search_vtx->axe = axe;
	search_vtx->alpha = alpha;
	search_vtx->N_half = N_half;
	search_vtx->using_bins = (100 < N);
	if (search_vtx->using_bins) {
		double bins_size =
		  (vertices[N-1]->x[axe] - vertices[0]->x[axe])/ sqrt(N);
		search_vtx->bins = nb_bins2D_create(bins_size);
		for (uint32_t i = 0; i < N; i++)
			nb_bins2D_insert(search_vtx->bins, vertices[i]);    
	} else {
		search_vtx->bins = NULL;
	}
}

static bool set_first_trg_into_AFL(nb_tessellator2D__t *mesh,
				   search_vtx_t *search_vtx,
				   nb_container_t *AFL)
{
	bool trg_created = false;
	if (0 < search_vtx->N_half) {
		msh_trg_t* first_trg =	create_1st_trg(mesh, search_vtx);
		if (NULL != first_trg) {
			mesh_add_triangle(mesh, first_trg);
			nb_container_insert(AFL, first_trg->s1);
			nb_container_insert(AFL, first_trg->s2);
			nb_container_insert(AFL, first_trg->s3);		
			mesh->do_after_insert_trg(mesh);
			trg_created = true;
		} /* else [the points are collinear] */
	} /* else [if all points are collinear]*/
	return trg_created;
}

static void clear_search_vtx(search_vtx_t *search_vtx)
{
	if (search_vtx->using_bins)
		nb_bins2D_destroy(search_vtx->bins);
}

static uint32_t triangulate_wall(nb_tessellator2D__t *mesh,
				 const search_vtx_t *search_vtx,
				 nb_container_t *AFL_alpha,
				 nb_container_t *AFL_1,
				 nb_container_t *AFL_2)
{
	uint32_t n_trg_alpha = 0;
	while (nb_container_is_not_empty(AFL_alpha)) {
		msh_edge_t* edge = nb_container_delete_first(AFL_alpha);
		msh_trg_t* trg = create_trg(mesh, search_vtx, edge);
		if (NULL != trg) {
			mesh_add_triangle(mesh, trg);
			update_AFLs(trg, edge, AFL_alpha, AFL_1, AFL_2,
				    search_vtx);
			n_trg_alpha += 1;
			mesh->do_after_insert_trg(mesh);
		}
	}
	return n_trg_alpha;
}

static void update_AFLs(const msh_trg_t *const trg,
			const msh_edge_t *const edge,
			nb_container_t *AFL_alpha,
			nb_container_t *AFL_1,
			nb_container_t *AFL_2,
			const search_vtx_t *search_vtx)
{
	msh_edge_t* complement[2];
	mtrg_get_complement_edges(trg, edge, complement);
	for (int8_t s = 0; s < 2; s++) {
		nb_container_t *side_AFL = 
			get_side_AFL(complement[s], search_vtx, 
				     AFL_1, AFL_alpha, AFL_2);
		update_AFL(side_AFL, complement[s]);
	}

}

static uint32_t dewall(nb_tessellator2D__t* mesh)
{
	nb_container_type cnt_type = NB_HASH;
	uint32_t vtx_size = mesh->N_input_vtx * sizeof(msh_vtx_t*);
	uint32_t memsize = vtx_size +
		nb_container_get_memsize(cnt_type);
	char *memblock = nb_soft_allocate_mem(memsize);

	msh_vtx_t **vertices = (void*) memblock;
	memcpy(vertices, mesh->input_vtx, vtx_size);

	nb_container_t* AFL = (void*) (memblock + vtx_size);
	nb_container_init(AFL, cnt_type);

	nb_container_set_key_generator(AFL, hash_key_edge);
	uint32_t N_trg = dewall_recursion(mesh, mesh->N_input_vtx,
					  vertices, 0, AFL);
  	nb_container_finish(AFL);
	nb_soft_free_mem(memsize, memblock);
	return N_trg;
}
