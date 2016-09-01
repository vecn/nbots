#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <alloca.h>

#include "nb/memory_bot.h"
#include "nb/container_bot.h"
#include "nb/geometric_bot/utils2D.h"
#include "nb/geometric_bot/knn/bins2D.h"
#include "nb/geometric_bot/knn/bins2D_iterator.h"
#include "nb/geometric_bot/mesh/mesh2D.h"

#include "nb/geometric_bot/mesh/modules2D/area_analizer.h"

#include "../mesh2D_structs.h"

#define POW2(a) ((a)*(a))

typedef struct {
	double area;
	nb_container_t *trgs;
} subarea_t;

static void* subarea_create(void);
static void subarea_destroy(void *subarea_ptr);
static void subarea_clear(void *subarea_ptr);
static int8_t subarea_compare_size(const void *subarea1_ptr,
				   const void *subarea2_ptr);

static void group_trg_by_areas(const nb_mesh_t *const mesh,
			       nb_container_t *areas);
static double* get_centroids(const nb_mesh_t *mesh,
			     nb_container_t *areas, uint32_t *N_centroids);
static double* get_centroids_if_enclosed(const nb_mesh_t *mesh,
				      nb_container_t *areas,
				      uint32_t *N_centroids);
static bool area_is_enclosed(const nb_container_t *area_trg);
static void calculate_area_centroid(const nb_mesh_t *mesh,
				    nb_container_t *area_trg,
				    double centroid[2]);
static msh_edge_t* check_if_internal_edge_is_longer(msh_edge_t *edge,
						    msh_edge_t *longest_edge,
						    double *max_length2);
static double spread_infection(msh_trg_t* trg_infected,
	       /* NULL if not required */ nb_container_t* infected_trg,
					  bool blocking_with_input_segments);

static int8_t compare_area1_isGreaterThan_area2(const void *const  a1,
						const void *const  a2);
static void get_useful_vtx(const nb_mesh_t *mesh, nb_container_t *useful_vtx);
static void delete_unused_vtx(nb_mesh_t *mesh, nb_container_t *useful_vtx);
static void get_unused_vtx(const nb_mesh_t *mesh,
			   const nb_container_t *useful_vtx,
			   nb_container_t *unused_vtx);
static void update_input_array(nb_mesh_t *mesh, const msh_vtx_t *vtx);
static uint16_t get_N_areas(const nb_mesh_t *mesh,
			    bool block_with_input_sgm);
static uint16_t count_areas_by_infection(nb_mesh_t *mesh,
					 bool block_with_input_sgm);
static void uninfect(nb_mesh_t *mesh);
static void get_area_ids(const nb_mesh_t *mesh, nb_container_t *areas,
			 uint16_t *trg_area_id);
static void set_id_to_trg_in_area(const nb_mesh_t *mesh,
				  nb_container_t *area_trg, uint16_t area_id,
				  uint16_t *trg_area_id);

static void* subarea_create(void)
{
	uint16_t size = sizeof(subarea_t) +
		nb_container_get_memsize(NB_QUEUE);
	char *memblock = malloc(size);
	subarea_t *subarea = (void*) memblock;
	subarea->trgs = (void*) (memblock + sizeof(subarea_t));
	nb_container_init(subarea->trgs, NB_QUEUE);
	return subarea;
}

static void subarea_destroy(void *subarea_ptr)
{
	subarea_t *subarea = subarea_ptr;
	nb_container_finish(subarea->trgs);
	free(subarea);
}

static void subarea_clear(void *subarea_ptr)
{
	subarea_t *subarea = subarea_ptr;
	while (nb_container_is_not_empty(subarea->trgs)) {
		msh_trg_t *trg = nb_container_delete_first(subarea->trgs);
		trg->status = CLEAN;
	}
}

static int8_t subarea_compare_size(const void *subarea1_ptr,
				   const void *subarea2_ptr)
{
	const subarea_t *subarea1 = subarea1_ptr;
	const subarea_t *subarea2 = subarea2_ptr;
	int8_t out;
	if (subarea1->area - subarea2->area < -NB_GEOMETRIC_TOL)
		out = -1;
	else if (subarea1->area - subarea2->area > NB_GEOMETRIC_TOL)
		out = 1;
	else
		out = 0;
	return out;
}

double* vcn_mesh_get_centroids_of_subareas(const nb_mesh_t *const mesh,
					   uint32_t* N_centroids)
{
	nb_container_t* areas = alloca(nb_container_get_memsize(NB_SORTED));
	nb_container_init(areas, NB_SORTED);
	nb_container_set_comparer(areas, subarea_compare_size);

	group_trg_by_areas(mesh, areas);

	double *centroids = get_centroids(mesh, areas, N_centroids);

	nb_container_finish(areas);

	return centroids;
}

static void group_trg_by_areas(const nb_mesh_t *const mesh,
			       nb_container_t *areas)
{
	nb_iterator_t* iter = alloca(nb_iterator_get_memsize());
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, mesh->ht_trg);
	while (nb_iterator_has_more(iter)) {
		msh_trg_t* trg = (msh_trg_t*) nb_iterator_get_next(iter);
		if (CLEAN == trg->status) {
			subarea_t *subarea = subarea_create();
			subarea->area =
				spread_infection(trg, subarea->trgs, true);
			nb_container_insert(areas, subarea);
		}
	}
	nb_iterator_finish(iter);

	uninfect((nb_mesh_t*)mesh);
}

static double* get_centroids(const nb_mesh_t *mesh,
			     nb_container_t *areas, uint32_t *N_centroids)
{
	*N_centroids = nb_container_get_length(areas);
	double *centroids = NULL;
	if (0 < *N_centroids) {
		centroids = calloc(*N_centroids * 2, sizeof(*centroids));
		uint32_t i = 0;
		while (nb_container_is_not_empty(areas)) {
			subarea_t *subarea = nb_container_delete_first(areas);
			calculate_area_centroid(mesh, subarea->trgs,
						&(centroids[i*2]));
			subarea_destroy(subarea);
			i += 1;
		}
	}
	return centroids;
}


static void calculate_area_centroid(const nb_mesh_t *mesh,
				    nb_container_t *area_trg,
				    double centroid[2])
{
	uint32_t N_trg = nb_container_get_length(area_trg);
	if (0 == N_trg) {
		memset(centroid, 0, 2 * sizeof(*centroid));
	} else if (1 == N_trg) {
		msh_trg_t* trg = nb_container_delete_first(area_trg);
		vcn_utils2D_trg_get_centroid(trg->v1->x, trg->v2->x,
					     trg->v3->x, centroid);
	} else {
		msh_edge_t *max_edge = NULL;
		double max_length2 = 0;
		while (nb_container_is_not_empty(area_trg)) {
			msh_trg_t* trg = nb_container_delete_first(area_trg);
	
			max_edge = check_if_internal_edge_is_longer(trg->s1, max_edge,
								    &max_length2);
			max_edge = check_if_internal_edge_is_longer(trg->s2, max_edge,
								    &max_length2);
			max_edge = check_if_internal_edge_is_longer(trg->s3, max_edge,
								    &max_length2);
		}
		centroid[0] = (max_edge->v1->x[0] + max_edge->v2->x[0]) / 2.0;
		centroid[1] = (max_edge->v1->x[1] + max_edge->v2->x[1]) / 2.0;
	}
	centroid[0] = centroid[0] / mesh->scale + mesh->xdisp;
	centroid[1] = centroid[1] / mesh->scale + mesh->ydisp;
}

static msh_edge_t* check_if_internal_edge_is_longer(msh_edge_t *edge,
						    msh_edge_t *longest_edge,
						    double *max_length2)
{
	msh_edge_t *longest = longest_edge;
	if (!medge_is_subsgm(edge)) {
		double dist2 =
			vcn_utils2D_get_dist2(edge->v1->x, edge->v2->x);
		if (dist2 > *max_length2) {
			*max_length2 = dist2;
			longest = edge;
		}
	}
	return longest;
}

static double spread_infection(msh_trg_t* trg_infected,
    /* NULL if not required */ nb_container_t* infected_trg,
			       bool blocking_with_input_segments)
{
	double area = 0.0;
	if (CLEAN == trg_infected->status) {
		trg_infected->status = INFECTED;

		if (NULL != infected_trg)
			nb_container_insert(infected_trg, trg_infected);

		bool segment_nonblock = !blocking_with_input_segments;
		if (blocking_with_input_segments)
			segment_nonblock = !medge_is_subsgm(trg_infected->s1);

		if (trg_infected->t1 != NULL && segment_nonblock) {
			msh_trg_t* nb_trg = trg_infected->t1;
			area += spread_infection(nb_trg,
						 infected_trg,
						 blocking_with_input_segments);
		}

		if (blocking_with_input_segments)
			segment_nonblock = !medge_is_subsgm(trg_infected->s2);

		if (NULL != trg_infected->t2 && segment_nonblock) {
			msh_trg_t* nb_trg = trg_infected->t2;
			area += spread_infection(nb_trg,
						 infected_trg,
						 blocking_with_input_segments);
		}

		if (blocking_with_input_segments)
			segment_nonblock = !medge_is_subsgm(trg_infected->s3);

		if (NULL != trg_infected->t3 && segment_nonblock) {
			msh_trg_t* nb_trg = trg_infected->t3;
			area += spread_infection(nb_trg,
						 infected_trg,
						 blocking_with_input_segments);
		}

		area += vcn_utils2D_get_trg_area(trg_infected->v1->x,
						 trg_infected->v2->x,
						 trg_infected->v3->x);
	}
	return area;
}

static int8_t compare_area1_isGreaterThan_area2
                                (const void *const restrict a1,
				 const void *const restrict a2)
{
	double area1_d = ((double*)((void**)a1)[0])[0];
	double area2_d = ((double*)((void**)a2)[0])[0];
	uint64_t area1 = (uint64_t)(1e8 * area1_d);
	uint64_t area2 = (uint64_t)(1e8 * area2_d);
	if (area2 < area1) 
		return 1;
	else if (area2 > area1)
		return -1;
	else
		return 0;
}

double* vcn_mesh_get_centroids_of_enveloped_areas(const nb_mesh_t *const mesh,
						  uint32_t* N_centroids)
{
	nb_container_t* areas = alloca(nb_container_get_memsize(NB_SORTED));
	nb_container_init(areas, NB_SORTED);
	nb_container_set_comparer(areas, subarea_compare_size);

	group_trg_by_areas(mesh, areas);

	double *centroids =
		get_centroids_if_enclosed(mesh, areas, N_centroids);

	nb_container_finish(areas);

	return centroids;
}

static double* get_centroids_if_enclosed(const nb_mesh_t *mesh,
					 nb_container_t *areas,
					 uint32_t *N_centroids)
{
	*N_centroids = nb_container_get_length(areas);
	
	double *out = NULL;
	if (0 < *N_centroids) {
		const uint32_t memsize = *N_centroids * 2 * sizeof(double);
		double *centroids = NB_SOFT_MALLOC(memsize);

		uint32_t i = 0;
		while (nb_container_is_not_empty(areas)) {
			subarea_t *subarea = nb_container_delete_first(areas);
			if (area_is_enclosed(subarea->trgs)) {
				calculate_area_centroid(mesh, subarea->trgs,
							&(centroids[i*2]));
				i += 1;
			} else {
				subarea_clear(subarea);
			}
			subarea_destroy(subarea);
		}
		*N_centroids = i;
		if (0 < *N_centroids) {
			uint32_t imemsize = i * 2 * sizeof(double);
			out = malloc(imemsize);
			memcpy(out, centroids, imemsize);
		}
		NB_SOFT_FREE(memsize, centroids);
	}
	return out;

}

static bool area_is_enclosed(const nb_container_t *area_trg)
{
	nb_iterator_t *iter = alloca(nb_iterator_get_memsize());
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, area_trg);
	bool out = true;
	while (nb_iterator_has_more(iter)) {
		const msh_trg_t *trg = nb_iterator_get_next(iter);
		if (NULL == trg->t1 || NULL == trg->t2 || NULL == trg->t3) {
			out = false;
			break;
		}
	}
	nb_iterator_finish(iter);
	return out;
}

double vcn_mesh_clear_enveloped_areas(nb_mesh_t* mesh,
				     double* area_removed)
{
	nb_container_t* areas = alloca(nb_container_get_memsize(NB_SORTED));
	nb_container_init(areas, NB_SORTED);

	nb_container_set_comparer(areas, compare_area1_isGreaterThan_area2);
	nb_iterator_t* iter = alloca(nb_iterator_get_memsize());
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, mesh->ht_trg);
	while (nb_iterator_has_more(iter)) {
		msh_trg_t* trg = (msh_trg_t*) nb_iterator_get_next(iter);
		if (CLEAN == trg->status) {
			nb_container_t* area_trg =
				nb_container_create(NB_SORTED);
			double* area = malloc(sizeof(*area));
			area[0] = spread_infection(trg, area_trg, true);
			void** obj = malloc(2 * sizeof(*obj));
			obj[0] = area;
			obj[1] = area_trg;
			nb_container_insert(areas, obj);
		}
	}
	nb_iterator_finish(iter);

	/* Destroy enveloped areas */
	if (NULL != area_removed)
		area_removed[0] = 0;
	while (nb_container_get_length(areas) > 1) {
		void** obj = nb_container_delete_first(areas);
		if(NULL != area_removed)
			area_removed[0] += ((double*)obj[0])[0];
		free(obj[0]);
		nb_container_t* area_trg = obj[1];
		while (nb_container_is_not_empty(area_trg)) {
			msh_trg_t* trg = nb_container_delete_first(area_trg);
			nb_container_delete(mesh->ht_trg, trg);
			mesh_substract_triangle(mesh, trg);
			mtrg_free(mesh, trg);
		}
		nb_container_destroy(area_trg);
		free(obj);
	}

	void** main_obj = nb_container_delete_first(areas);
	double ret_val = ((double*)main_obj[0])[0];
	free(main_obj[0]);
	nb_container_t* main_area = (nb_container_t*)main_obj[1];
	while (nb_container_is_not_empty(main_area)) {
		msh_trg_t* trg = nb_container_delete_first(main_area);
		trg->status = CLEAN;
	}
	nb_container_destroy(main_area);
	free(main_obj);
	/* Free memory */
	nb_container_finish(areas);
  
	/* Rescale value */
	ret_val /= POW2(mesh->scale);
	if(NULL != area_removed)
		area_removed[0] /= POW2(mesh->scale);

	/* Return biggest area */
	return ret_val;
}

double vcn_mesh_keep_biggest_continuum_area(nb_mesh_t* mesh,
					    double* area_removed)
{
	nb_container_t* areas = alloca(nb_container_get_memsize(NB_SORTED));
	nb_container_init(areas, NB_SORTED);

	nb_container_set_comparer(areas, compare_area1_isGreaterThan_area2);
	nb_iterator_t* iter = alloca(nb_iterator_get_memsize());
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, mesh->ht_trg);
	while (nb_iterator_has_more(iter)) {
		msh_trg_t* trg = (msh_trg_t*)nb_iterator_get_next(iter);
		if (CLEAN == trg->status) {
			nb_container_t* area_trg =
				nb_container_create(NB_SORTED);
			double* area = malloc(sizeof(*area));
			area[0] = spread_infection(trg, area_trg, false);
			void** obj = malloc(2*sizeof(*obj));
			obj[0] = area;
			obj[1] = area_trg;
			nb_container_insert(areas, obj);
		}
	}
	nb_iterator_finish(iter);

	/* Remove separated areas */
	if (NULL != area_removed)
		area_removed[0] = 0;
	while (nb_container_get_length(areas) > 1) {
		void** obj = nb_container_delete_first(areas);
		if (NULL != area_removed) 
			area_removed[0] += ((double*)obj[0])[0];
		free(obj[0]);
		nb_container_t* area_trg = (nb_container_t*)obj[1];
		while (nb_container_is_not_empty(area_trg)) {
			msh_trg_t* trg = nb_container_delete_first(area_trg);
			nb_container_delete(mesh->ht_trg, trg);
			mesh_substract_triangle(mesh, trg);
			mtrg_free(mesh, trg);
		}
		nb_container_destroy(area_trg);
		free(obj);
	}
	/* Keep biggest area */
	void** main_obj = nb_container_delete_first(areas);
	double ret_val = ((double*)main_obj[0])[0];
	free(main_obj[0]);
	nb_container_t* main_area = (nb_container_t*)main_obj[1];
	while (nb_container_is_not_empty(main_area)) {
		msh_trg_t* trg = nb_container_delete_first(main_area);
		trg->status = CLEAN;
	}
	nb_container_destroy(main_area);
	free(main_obj);
	nb_container_finish(areas);

	/* Rescale value */
	ret_val /= POW2(mesh->scale);
	if (NULL != area_removed)
		area_removed[0] /= POW2(mesh->scale);

	/* Return biggest area */
	return ret_val;
}

uint32_t vcn_mesh_delete_isolated_segments(nb_mesh_t *const restrict mesh)
{
	uint32_t removed = 0;
	for (uint32_t i = 0; i < mesh->N_input_sgm; i++) {
		msh_edge_t* restrict sgm = mesh->input_sgm[i];
		if (NULL != sgm) {
			if (NULL == sgm->t1 && NULL == sgm->t2) {
				mesh->input_sgm[i] = NULL;
				while (NULL != sgm) {
					msh_edge_t* to_free = sgm;
					sgm = medge_subsgm_next(sgm);
					nb_container_delete(mesh->ht_edge,
							    to_free);
					medge_destroy_subsgm_attribute(to_free);
					medge_free(mesh, to_free);
					removed += 1;
				}
			}
		}
	}
	return removed;
}

uint32_t vcn_mesh_delete_internal_input_segments(nb_mesh_t *const restrict mesh)
{
	uint32_t removed = 0;
	for (uint32_t i = 0; i < mesh->N_input_sgm; i++) {
		msh_edge_t* restrict sgm = mesh->input_sgm[i];
		if (NULL != sgm) {
			if (NULL != sgm->t1 && NULL != sgm->t2) {
				mesh->input_sgm[i] = NULL;      
				while (NULL != sgm) {
					msh_edge_t* prev_sgm = sgm;
					sgm = medge_subsgm_next(sgm);
					medge_destroy_subsgm_attribute(prev_sgm);
					removed += 1;
				}
			}
		}
	}
	return removed;
}

uint32_t vcn_mesh_delete_isolated_vertices(nb_mesh_t* mesh)
{
	nb_container_t* useful_vtx =
		alloca(nb_container_get_memsize(NB_SORTED));
	nb_container_init(useful_vtx, NB_SORTED);

	get_useful_vtx(mesh, useful_vtx);	
	uint32_t length = vcn_bins2D_get_length(mesh->ug_vtx);
	delete_unused_vtx(mesh, useful_vtx);

	nb_container_finish(useful_vtx);

	return length - vcn_bins2D_get_length(mesh->ug_vtx);
}

static void get_useful_vtx(const nb_mesh_t *mesh, nb_container_t *useful_vtx)
{
	nb_iterator_t* iter = alloca(nb_iterator_get_memsize());
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, mesh->ht_trg);
	while (nb_iterator_has_more(iter)) {
		msh_trg_t* trg = (msh_trg_t*)nb_iterator_get_next(iter);
		nb_container_insert(useful_vtx, trg->v1);
		nb_container_insert(useful_vtx, trg->v2);
		nb_container_insert(useful_vtx, trg->v3);
	}
	nb_iterator_finish(iter);
}

static void delete_unused_vtx(nb_mesh_t *mesh, nb_container_t *useful_vtx)
{
	nb_container_t *unused_vtx = alloca(nb_container_get_memsize(NB_QUEUE));
	nb_container_init(unused_vtx, NB_QUEUE);

	get_unused_vtx(mesh, useful_vtx, unused_vtx);

	while (nb_container_is_not_empty(unused_vtx)) {
		msh_vtx_t *vtx = nb_container_delete_first(unused_vtx);
		vcn_bins2D_delete(mesh->ug_vtx, vtx);
		update_input_array(mesh, vtx);
		mvtx_destroy(mesh, vtx);
	}
	nb_container_finish(unused_vtx);
}

static void get_unused_vtx(const nb_mesh_t *mesh,
			   const nb_container_t *useful_vtx,
			   nb_container_t *unused_vtx)
{
	vcn_bins2D_iter_t* iter = alloca(vcn_bins2D_iter_get_memsize());
	vcn_bins2D_iter_init(iter);
	vcn_bins2D_iter_set_bins(iter, mesh->ug_vtx);
	while (vcn_bins2D_iter_has_more(iter)) {
		const msh_vtx_t* vtx = vcn_bins2D_iter_get_next(iter);
		if (!nb_container_exist(useful_vtx, vtx)) {
			nb_container_insert(unused_vtx, vtx);
		}
	}
	vcn_bins2D_iter_finish(iter);
}

static void update_input_array(nb_mesh_t *mesh, const msh_vtx_t *vtx)
{
	if (mvtx_is_type_origin(vtx, INPUT)) {
		for (uint32_t i = 0; i < mesh->N_input_vtx; i++) {
			if (vtx == mesh->input_vtx[i]) {
				mesh->input_vtx[i] = NULL;
				goto EXIT;
			}
		}
	}
EXIT:
	return;
}

bool vcn_mesh_is_continuum(const nb_mesh_t *mesh)
{
	uint16_t N = vcn_mesh_get_N_continuum_areas(mesh);
	return (N == 1);
}

inline uint16_t vcn_mesh_get_N_subareas(const nb_mesh_t *mesh)
{
	return get_N_areas(mesh, true);
}

static uint16_t get_N_areas(const nb_mesh_t *mesh,
			    bool block_with_input_sgm)
{
	uint16_t counter = count_areas_by_infection((nb_mesh_t*)mesh,
						    block_with_input_sgm);
	uninfect((nb_mesh_t*)mesh);
	return counter;
}

static uint16_t count_areas_by_infection(nb_mesh_t *mesh,
					 bool block_with_input_sgm)
{
	nb_iterator_t* iter = alloca(nb_iterator_get_memsize());
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, mesh->ht_trg);

	uint16_t counter = 0;
	while (nb_iterator_has_more(iter)) {
		msh_trg_t* trg = (msh_trg_t*) nb_iterator_get_next(iter);
		if (CLEAN == trg->status) {
			spread_infection(trg, NULL, block_with_input_sgm);
			counter += 1;
		}
	}
	nb_iterator_finish(iter);
	return counter;
}

static void uninfect(nb_mesh_t *mesh)
{
	nb_iterator_t* iter = alloca(nb_iterator_get_memsize());
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, mesh->ht_trg);

	while (nb_iterator_has_more(iter)) {
		msh_trg_t* trg = (msh_trg_t*) nb_iterator_get_next(iter);
		trg->status = CLEAN;
	}
	nb_iterator_finish(iter);	
}

uint16_t nb_mesh_get_subareas(const nb_mesh_t *mesh, uint16_t *area_id)
{	
	nb_container_t* areas = alloca(nb_container_get_memsize(NB_SORTED));
	nb_container_init(areas, NB_SORTED);
	nb_container_set_comparer(areas, subarea_compare_size);

	group_trg_by_areas(mesh, areas);
	
	uint16_t N_areas = nb_container_get_length(areas);

	get_area_ids(mesh, areas, area_id);

	nb_container_finish(areas);

	return N_areas;
}

static void get_area_ids(const nb_mesh_t *mesh, nb_container_t *areas,
			 uint16_t *trg_area_id)
{
	uint16_t area_id = 0;
	while (nb_container_is_not_empty(areas)) {
		subarea_t *subarea = nb_container_delete_first(areas);
		set_id_to_trg_in_area(mesh, subarea->trgs, area_id, trg_area_id);
		subarea_destroy(subarea);
		area_id += 1;
	}
}

static void set_id_to_trg_in_area(const nb_mesh_t *mesh,
				  nb_container_t *area_trg, uint16_t area_id,
				  uint16_t *trg_area_id)
{
	while (nb_container_is_not_empty(area_trg)) {
		msh_trg_t* trg = nb_container_delete_first(area_trg);
		uint32_t tid = trg->id;
		trg_area_id[tid] = area_id;
	}

}

inline uint16_t vcn_mesh_get_N_continuum_areas(const nb_mesh_t *mesh)
{
	return get_N_areas(mesh, false);
}
