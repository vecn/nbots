#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include "nb/geometric_bot/point2D.h"
#include "nb/geometric_bot/utils2D.h"
#include "nb/geometric_bot/knn/bins2D.h"
#include "nb/geometric_bot/mesh/tessellator2D.h"
#include "nb/geometric_bot/mesh/dewall.h"

#include "tessellator2D_structs.h" // TODO: Remove dependency here

#include "dewall_dependencies.h"

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
	vtx_t** vtx_array;
	int8_t axe;
	double alpha;
	int N_half;
} search_vtx_t;

static  void mesh_init_disp_and_scale
                            (nb_tessellator2D_t *const mesh,
			     uint32_t N_vertices,
			     const double *const vertices);
static int8_t compare_using_axe(const void *const v1_ptr,
				const void *const v2_ptr,
				const void *const axe_ptr);
static vtx_t* find_1st_vtx(search_vtx_t *search_vtx);
static vtx_t* find_2nd_vtx(search_vtx_t* search_vtx, vtx_t* v1);
static vtx_t* find_2nd_vtx_exahustive_search(uint32_t N, 
						vtx_t** vertices,
						vtx_t* v1);
static bool compare_2nd_vtx(const nb_point2D_t *const p_ref,
			    const nb_point2D_t *const p,
			    const void *const data);
static vtx_t* find_2nd_vtx_using_bins(nb_bins2D_t* bins, vtx_t* v1,
					 int8_t axe, double alpha);
static double swap_vtx_if_minimize_delaunay_distance
                                      (const mesh_t *const mesh,
				       const vtx_t *const v1,
				       const vtx_t *const v2,
				       vtx_t** v3, double min_dist,
				       vtx_t* vtx, double dist);
static bool is_cocircular(const mesh_t *const mesh,
			  const vtx_t *const v1,
			  const vtx_t *const v2,
			  const vtx_t *const v3);
static vtx_t* find_3rd_vtx(const mesh_t *const mesh,
			      const search_vtx_t *const search_vtx,
			      const vtx_t *const  v1,
			      const vtx_t *const  v2);
static vtx_t* find_3rd_vtx_exahustive_search(const mesh_t *const mesh,
						uint32_t N, 
						vtx_t** vertices,
						const vtx_t *const v1,
						const vtx_t *const v2);
static double set_v3(const mesh_t *const mesh,
		     const vtx_t *const v1, const vtx_t *const v2,
		     vtx_t* v3_candidate, vtx_t** v3, double min_dist);
static void init_bins2D_vertices_param(nb_bins2D_vertices_t *vtx_queue,
				       queue_t *vertices);
static vtx_t* find_3rd_vtx_using_bins(const mesh_t *const mesh,
					 const nb_bins2D_t *const bins,
					 const vtx_t *const  v1,
					 const vtx_t *const  v2);
static trg_t* create_1st_trg(mesh_t *mesh,
				 search_vtx_t* search_vtx);
static trg_t* create_trg(mesh_t *mesh,
			     const search_vtx_t *search_vtx,
			     msh_edge_t *edge);
static void update_AFL(afl_t *AFL, const msh_edge_t *const edge);
static uint32_t dewall_recursion
                       (mesh_t *mesh, uint32_t N,
			vtx_t** vertices,
			uint16_t deep_level,
			afl_t* AFL);
static uint32_t triangulate_wall(mesh_t *mesh,
				 const search_vtx_t *search_vtx,
				 afl_t *AFL_alpha,
				 afl_t *AFL_1,
				 afl_t *AFL_2);
static void update_AFLs(const trg_t *const trg,
			const msh_edge_t *const edge,
			afl_t *AFL_alpha,
			afl_t *AFL_1,
			afl_t *AFL_2,
			const search_vtx_t *search_vtx);
static double get_alpha(uint32_t N, vtx_t **vertices,
			int8_t axe, uint32_t N_mid);
static afl_t* select_side_AFL(const msh_edge_t *const edge,
				     const search_vtx_t *search_vtx,
				     afl_t *AFL_neg,
				     afl_t *AFL_zero,
				     afl_t *AFL_pos);

static uint32_t split_vtx_array(const afl_t *const AFL,
				uint32_t N, vtx_t **vertices,
				int8_t axe);
static void set_interval(interval_t *interval,
			 const afl_t *const AFL,
			 uint32_t N, vtx_t **vertices,
			 int8_t axe);
static uint32_t get_mid_inside_interval(interval_t *interval, uint32_t N,
					vtx_t **vertices, int8_t axe);
static bool not_space_for_alpha(vtx_t **vertices, uint32_t i, int8_t axe);
static void init_search_vtx(search_vtx_t *search_vtx, uint32_t N,
			    vtx_t **vertices, int8_t axe,
			    double alpha, uint32_t N_half);
static bool set_first_trg_into_AFL(mesh_t *mesh,
				   search_vtx_t *search_vtx, afl_t *AFL);
static void clear_search_vtx(search_vtx_t *search_vtx);

static uint32_t dewall_memsize(uint32_t n_vtx);
static uint32_t dewall(mesh_t* mesh, uint32_t n_vtx, char *memblock);

void nb_tessellator2D_get_delaunay(nb_tessellator2D_t *mesh, uint32_t N_vertices,
				   const double *const vertices)
{
	nb_tessellator2D_clear(mesh);
	if (0 == N_vertices)
		goto EXIT;

	mesh_init_disp_and_scale(mesh, N_vertices, vertices);
	mesh->N_input_vtx = N_vertices;
	mesh->input_vtx =
		module()->mem.allocate_zero(mesh->N_input_vtx *
					    sizeof(*(mesh->input_vtx)));

	for (uint32_t i = 0; i < N_vertices; i++) {
		vtx_t* vtx = mvtx_create(mesh);
		mesh->input_vtx[i] = vtx;
		mvtx_set_type_origin(vtx, INPUT);
		vtx->x[0] = mesh->scale *
			(vertices[i * 2] - mesh->xdisp);
		vtx->x[1] = mesh->scale *
			(vertices[i*2+1] - mesh->ydisp);
		nb_bins2D_insert(mesh->ug_vtx, vtx);
	}
	uint32_t memsize = dewall_memsize(N_vertices);
	char *memblock = module()->mem.allocate(memsize);

	dewall((mesh_t*)mesh, N_vertices, memblock);

	module()->mem.free(memblock);
EXIT:
	return;
}

static inline void mesh_init_disp_and_scale
                            (nb_tessellator2D_t *const restrict mesh,
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

static int8_t compare_using_axe(const void *const v1_ptr,
				const void *const v2_ptr,
				const void *const axe_ptr)
{
	const vtx_t *restrict v1 = *((vtx_t**)v1_ptr);
	const vtx_t *restrict v2 = *((vtx_t**)v2_ptr);
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

static vtx_t* find_2nd_vtx(search_vtx_t* search_vtx, vtx_t* v1)
{
	vtx_t *v2;
	if (search_vtx->using_bins)
		v2 = find_2nd_vtx_using_bins(search_vtx->bins, v1,
					    search_vtx->axe,
					    search_vtx->alpha);
	else
		v2 = find_2nd_vtx_exahustive_search(search_vtx->N_half,
						   search_vtx->vtx_array, v1);
	return v2;
}

static vtx_t* find_2nd_vtx_exahustive_search
                                               (uint32_t N, 
						vtx_t** vertices,
						vtx_t* v1)
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

static bool compare_2nd_vtx(const nb_point2D_t *const restrict p_ref,
				   const nb_point2D_t *const restrict p,
				   const void *const restrict data)
{
	register double alpha = ((double*)data)[0];
	register int8_t axe = (uint8_t)((double*)data)[1];
	return (alpha - p->x[axe] > 0.0);

}

static vtx_t* find_2nd_vtx_using_bins(nb_bins2D_t* bins,
						vtx_t* v1, int8_t axe,
						double alpha)
{
	double filter_data[2];
	filter_data[0] = alpha;
	filter_data[1] = (double) axe;

	vtx_t* v2 = NULL;
	double dist;
	nb_bins2D_set_filter(bins, compare_2nd_vtx);
	nb_bins2D_set_filter_data(bins, filter_data);
	nb_bins2D_get_knn(bins, v1, 1, &v2, &dist);
	return v2;
}

static double swap_vtx_if_minimize_delaunay_distance
                                      (const mesh_t *const mesh,
				       const vtx_t *const v1,
				       const vtx_t *const v2,
				       vtx_t** v3, double min_dist,
				       vtx_t *vtx, double dist)
{
	bool do_swap = false;
	if (fabs(dist - min_dist) < NB_GEOMETRIC_TOL) {
	        do_swap = is_cocircular(mesh, v1, v2, *v3);
	} else {
		do_swap = dist < min_dist;
	}
	
	if (do_swap) {
		*v3 = vtx;
		min_dist = dist;
	}
	return min_dist;
}

static bool is_cocircular(const mesh_t *const mesh,
			  const vtx_t *const v1,
			  const vtx_t *const v2,
			  const vtx_t *const v3)
{
	return  module()->mesh.is_v3_intersecting_any_edge(mesh, v1, v2, v3);
}

static vtx_t* find_3rd_vtx(const mesh_t *const mesh,
				     const search_vtx_t *const search_vtx,
				     const vtx_t *const  v1,
				     const vtx_t *const  v2)
{
	vtx_t *v3;
	if (search_vtx->using_bins)
		v3 = find_3rd_vtx_using_bins(mesh, search_vtx->bins, v1, v2);
	else
		v3 = find_3rd_vtx_exahustive_search(mesh, search_vtx->N,
						   search_vtx->vtx_array,
						   v1, v2);
	return v3;
}

static vtx_t* find_3rd_vtx_exahustive_search(const mesh_t *const mesh,
						uint32_t N, 
						vtx_t** vertices,
						const vtx_t *const v1,
						const vtx_t *const v2)
{
	vtx_t *v3 = NULL;
	double min_dist = 0.0;
	for (uint32_t i = 0; i < N; i++) {
		if (vertices[i] != v1 && vertices[i] != v2) {
			if (nb_utils2D_is_in_half_side(v1->x, v2->x,
							vertices[i]->x))
				min_dist = set_v3(mesh, v1, v2, vertices[i],
						  &v3, min_dist);
		}
	}
	return v3;
}

static double set_v3(const mesh_t *const mesh,
		     const vtx_t *const v1, const vtx_t *const v2,
		     vtx_t* v3_candidate, vtx_t** v3, double min_dist)
{
	double dist = nb_utils2D_get_delaunay_dist(v1->x, v2->x,
						    v3_candidate->x);
	if (NULL == *v3) {
		/* Initialize minimum distance */
		min_dist = dist;
		*v3 = v3_candidate;
	} else {
		min_dist = 
			swap_vtx_if_minimize_delaunay_distance(mesh,
							       v1, v2, v3,
							       min_dist,
							       v3_candidate,
							       dist);
	}
	return min_dist;

}

static void init_bins2D_vertices_param(nb_bins2D_vertices_t *vtx_queue,
				       queue_t *vertices)
{
	vtx_queue->vtx = vertices;
	vtx_queue->add =
		(void (*)(void*, const void *const))module()->queue.add;
	vtx_queue->is_empty =
		(bool (*)(const void *const))module()->queue.is_empty;
}

static vtx_t* find_3rd_vtx_using_bins(const mesh_t *const mesh,
					 const nb_bins2D_t *const restrict bins,
					 const vtx_t *const restrict v1,
					 const vtx_t *const restrict v2)
{
	uint32_t cnt_size = module()->queue.size();
	queue_t* vertices = module()->mem.allocate(cnt_size);
	module()->queue.init(vertices);

	nb_bins2D_vertices_t vtx_queue;
	init_bins2D_vertices_param(&vtx_queue, vertices);
	nb_bins2D_get_candidate_points_to_min_delaunay(bins, v1, v2,
						       &vtx_queue);
	vtx_t *v3 = NULL;
	double min_dist = 0.0;
	while (!module()->queue.is_empty(vertices)) {
		vtx_t* vtx = module()->queue.poll(vertices);
		min_dist = set_v3(mesh, v1, v2, vtx, &v3, min_dist);
	}
	module()->queue.finish(vertices);
	module()->mem.free(vertices);
       	return v3;
}

static trg_t* create_1st_trg(mesh_t *mesh, search_vtx_t* search_vtx)
{
	vtx_t *v1 = find_1st_vtx(search_vtx);
	vtx_t *v2 = find_2nd_vtx(search_vtx, v1);
	vtx_t *v3 = find_3rd_vtx(mesh, search_vtx, v1, v2);
  
	if (NULL == v3) {
		/* Check in the other halfspace */
		vtx_t *aux_vtx = v1;
		v1 = v2;
		v2 = aux_vtx;
		v3 = find_3rd_vtx(mesh, search_vtx, v1, v2);
	}

	if (NULL == v3) {
		/* All points are collinear */
		return NULL;
	}
	return module()->mesh.new_triangle(mesh, v1, v2, v3);
}

static vtx_t* find_1st_vtx(search_vtx_t *search_vtx)
{
	/* Closest vertex to the wall */
	return search_vtx->vtx_array[search_vtx->N_half];
}

static trg_t* create_trg(mesh_t *mesh,
			     const search_vtx_t *search_vtx,
			     msh_edge_t *edge)
{
	/* Select the correct segment orientation */
	vtx_t *restrict v1;
	vtx_t *restrict v2;
	if (NULL == edge->t1) {
		v1 = edge->v1;
		v2 = edge->v2;
	} else {
		v1 = edge->v2;
		v2 = edge->v1;
	}

	vtx_t *restrict v3 = find_3rd_vtx(mesh, search_vtx, v1, v2);

	if(NULL == v3) {
		return NULL;
	}
	return module()->mesh.new_triangle(mesh, v1, v2, v3);
}

static void update_AFL(afl_t *AFL, const msh_edge_t *const edge)
{
	if (NULL == module()->afl.delete(AFL, edge)) {
		module()->afl.insert(AFL, edge);
	}
}

static uint32_t dewall_recursion
                       (mesh_t *mesh, uint32_t N,
			vtx_t** vertices,
			uint16_t deep_level,
			afl_t* AFL)
{
	uint32_t N_trg = 0;
	if (N < 3)
		goto EXIT;

	/* Ascending sorting of vertices depending on the axe divisor */
	int8_t axe = deep_level % 2; /* 0:X, 1:Y */

	nb_qsort_wd(vertices, N, sizeof(*vertices),
		     compare_using_axe, &axe);
	uint32_t N_mid = split_vtx_array(AFL, N, vertices, axe);
	double alpha = get_alpha(N, vertices, axe, N_mid);

	search_vtx_t search_vtx;
	init_search_vtx(&search_vtx, N, vertices, axe, alpha, N_mid);

	/* Initialize triangles counter */
	uint32_t n_trg_alpha = 0;

	/* Create first triangle */
	if (module()->afl.is_empty(AFL)) {
		if (set_first_trg_into_AFL(mesh, &search_vtx, AFL))
			n_trg_alpha += 1;
		else
			goto EXIT;
	}

	/* Initialize Action face lists */
	uint32_t AFL_size = module()->afl.size();
	uint32_t memsize = 3 * AFL_size;
	char *memblock = module()->mem.allocate(memsize);

	afl_t* AFL_alpha = (void*) memblock;
	module()->afl.init(AFL_alpha, hash_key_edge);

	afl_t* AFL_1 = (void*) (memblock + AFL_size);
	module()->afl.init(AFL_1, hash_key_edge);

	afl_t* AFL_2 = (void*) (memblock + 2 * AFL_size);
	module()->afl.init(AFL_2, hash_key_edge);

	/* Redistribute segments from the main AFL to the three AFLs */
	while (!module()->afl.is_empty(AFL)) {
		msh_edge_t* edge = module()->afl.delete_any(AFL);
		afl_t *side_AFL = select_side_AFL(edge, &search_vtx,
						  AFL_1, AFL_alpha,
						  AFL_2);
		module()->afl.insert(side_AFL, edge);
	}
	
	n_trg_alpha += triangulate_wall(mesh, &search_vtx, AFL_alpha,
					AFL_1, AFL_2);

	module()->afl.finish(AFL_alpha);

	clear_search_vtx(&search_vtx);

	/* Recursive triangulation */
	uint32_t n_trg_1 = 0;
	if (!module()->afl.is_empty(AFL_1))
		n_trg_1 += dewall_recursion(mesh, N_mid, vertices,
					    deep_level + 1, AFL_1);
	module()->afl.finish(AFL_1);

	uint32_t n_trg_2 = 0;
	if (!module()->afl.is_empty(AFL_2))
		n_trg_2 += dewall_recursion(mesh, N - N_mid, &(vertices[N_mid]),
					    deep_level + 1, AFL_2);
	module()->afl.finish(AFL_2);

	N_trg = n_trg_alpha + n_trg_1 + n_trg_2;

	module()->mem.free(memblock);
EXIT:
	return N_trg;
}

static double get_alpha(uint32_t N, vtx_t **vertices,
			       int8_t axe, uint32_t N_mid)
{
	if (N > 3) {
		if (0 == N_mid)
			return vertices[0]->x[axe];
		else
			return 0.6 * vertices[N_mid]->x[axe] +
				0.4 * vertices[N_mid-1]->x[axe];
	} else {
		double min = MIN(vertices[0]->x[axe], vertices[1]->x[axe]);
		min = MIN(min, vertices[2]->x[axe]);

		double max = MAX(vertices[0]->x[axe], vertices[1]->x[axe]);
		max = MAX(max, vertices[2]->x[axe]);

		return (min + max)/2.0;
	}
}

static afl_t* select_side_AFL(const msh_edge_t *const edge,
				     const search_vtx_t *search_vtx,
				     afl_t *AFL_neg,
				     afl_t *AFL_zero,
				     afl_t *AFL_pos)
{
	bool v1 = (edge->v1->x[search_vtx->axe] - search_vtx->alpha > 0.0);
	bool v2 = (edge->v2->x[search_vtx->axe] - search_vtx->alpha > 0.0);
	if (v1 != v2)
		return AFL_zero;
	else if (v1)
		return AFL_pos;
	else 
		return AFL_neg;
}

static uint32_t split_vtx_array(const afl_t *const AFL,
				uint32_t N, vtx_t **vertices,
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
			 const afl_t *const AFL,
			 uint32_t N, vtx_t **vertices,
			 int8_t axe)
{
	if (module()->afl.is_empty(AFL)) {
		interval->min = vertices[0]->x[axe];
		interval->max = vertices[N-1]->x[axe];
	} else {
		double range[2];
		module()->afl.get_range_with_faces(AFL, axe, range);
		interval->min = range[0];
		interval->max = range[1];
	}
}

static uint32_t get_mid_inside_interval(interval_t *interval, uint32_t N,
					vtx_t **vertices, int8_t axe)
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

static bool not_space_for_alpha(vtx_t **vertices, uint32_t i,
				       int8_t axe)
{
	return (vertices[i]->x[axe] - vertices[i-1]->x[axe])
		< NB_GEOMETRIC_TOL;
}

static void init_search_vtx(search_vtx_t *search_vtx, uint32_t N,
			    vtx_t **vertices, int8_t axe,
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

static bool set_first_trg_into_AFL(mesh_t *mesh,
				   search_vtx_t *search_vtx,
				   afl_t *AFL)
{
	bool trg_created = false;
	if (0 < search_vtx->N_half) {
		trg_t* first_trg =	create_1st_trg(mesh, search_vtx);
		if (NULL != first_trg) {
			module()->afl.insert(AFL, first_trg->s1);
			module()->afl.insert(AFL, first_trg->s2);
			module()->afl.insert(AFL, first_trg->s3);
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

static uint32_t triangulate_wall(mesh_t *mesh,
				 const search_vtx_t *search_vtx,
				 afl_t *AFL_alpha,
				 afl_t *AFL_1,
				 afl_t *AFL_2)
{
	uint32_t n_trg_alpha = 0;
	while (!module()->afl.is_empty(AFL_alpha)) {
		msh_edge_t* edge = module()->afl.delete_any(AFL_alpha);
		trg_t* trg = create_trg(mesh, search_vtx, edge);
		if (NULL != trg) {
			update_AFLs(trg, edge, AFL_alpha, AFL_1, AFL_2,
				    search_vtx);
			n_trg_alpha += 1;
		}
	}
	return n_trg_alpha;
}

static void update_AFLs(const trg_t *const trg,
			const msh_edge_t *const edge,
			afl_t *AFL_alpha,
			afl_t *AFL_1,
			afl_t *AFL_2,
			const search_vtx_t *search_vtx)
{
	msh_edge_t* complement[2];
	mtrg_get_complement_edges(trg, edge, complement);
	for (int8_t s = 0; s < 2; s++) {
		afl_t *side_AFL = 
			select_side_AFL(complement[s], search_vtx,
					AFL_1, AFL_alpha, AFL_2);
		update_AFL(side_AFL, complement[s]);
	}

}

static uint32_t dewall_memsize(uint32_t n_vtx)
{
	uint32_t vtx_size = n_vtx * sizeof(vtx_t*);
	uint32_t afl_size = module()->afl.size();
	return vtx_size + afl_size;
}

static uint32_t dewall(mesh_t* mesh, uint32_t n_vtx, char *memblock)
{
	uint32_t vtx_size = n_vtx * sizeof(vtx_t*);

	vtx_t **vertices = (void*) memblock;

	memcpy(vertices,
	       ((nb_tessellator2D_t*)mesh)->input_vtx, vtx_size);// TODO

	afl_t *AFL = (void*) (memblock + vtx_size);
	module()->afl.init(AFL, hash_key_edge);

	uint32_t n_trg = dewall_recursion((mesh_t*)mesh, n_vtx,
					  vertices, 0, AFL);

	module()->afl.finish(AFL);
	return n_trg;
}
