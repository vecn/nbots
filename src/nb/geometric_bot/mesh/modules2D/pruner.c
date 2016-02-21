#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <alloca.h>

#include "nb/container_bot.h"
#include "nb/geometric_bot/utils2D.h"
#include "nb/geometric_bot/knn/bins2D.h"
#include "nb/geometric_bot/knn/bins2D_iterator.h"
#include "nb/geometric_bot/mesh/mesh2D.h"

#include "../mesh2D_structs.h"

#define _NB_SUBSEGMENT_VTX ((void*)0x2) /* REFACTOR: Taken from ruppert.c */

#define POW2(a) ((a)*(a))

static void group_trg_by_areas(const vcn_mesh_t *const mesh,
			       nb_container_t *areas);
static double* get_centroids(const vcn_mesh_t *mesh,
			     nb_container_t *areas, uint32_t *N_centroids);
static void calculate_area_centroid(const vcn_mesh_t *mesh,
				    nb_container_t *area_trg,
				    double centroid[2]);
static msh_edge_t* check_if_internal_edge_is_longer(msh_edge_t *edge,
						    msh_edge_t *longest_edge,
						    double *max_length2);
static double colorize_infected_triangles(msh_trg_t* trg_infected,
					  nb_container_t* l_infected_trg,
					  bool blocking_with_input_segments);

static int8_t compare_area1_isGreaterThan_area2(const void *const  a1,
						const void *const  a2);

double* vcn_mesh_get_centroids_of_subareas(const vcn_mesh_t *const mesh,
					   uint32_t* N_centroids)
{
	nb_container_t* areas = alloca(nb_container_get_memsize(NB_SORTED));
	nb_container_init(areas, NB_SORTED);
	nb_container_set_comparer(areas, compare_area1_isGreaterThan_area2);

	group_trg_by_areas(mesh, areas);

	double *centroids = get_centroids(mesh, areas, N_centroids);

	nb_container_finish(areas);

	return centroids;
}

static void group_trg_by_areas(const vcn_mesh_t *const mesh,
			       nb_container_t *areas)
{
	nb_iterator_t* iter = alloca(nb_iterator_get_memsize());
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, mesh->ht_trg);
	while (nb_iterator_has_more(iter)) {
		msh_trg_t* trg = (msh_trg_t*) nb_iterator_get_next(iter);
		if (NULL != trg->attr) {
			attr_t* attr = trg->attr;
			if(attr->id == 0x19FEC7ED)
				continue;
		}
		nb_container_t* area_trg = nb_container_create(NB_SORTED);
		double* area = malloc(sizeof(*area));
		*area = colorize_infected_triangles(trg, area_trg, true);
		void** obj = malloc(2 * sizeof(*obj));
		obj[0] = area;
		obj[1] = area_trg;
		nb_container_insert(areas, obj);
	}
	nb_iterator_finish(iter);
}

static double* get_centroids(const vcn_mesh_t *mesh,
			     nb_container_t *areas, uint32_t *N_centroids)
{
	*N_centroids = nb_container_get_length(areas);
	double *centroids = calloc(*N_centroids * 2, sizeof(*centroids));
	uint32_t i = 0;
	while (nb_container_is_not_empty(areas)) {
		void** obj = nb_container_delete_first(areas);
		nb_container_t *area_trg = obj[1];
		calculate_area_centroid(mesh, area_trg, &(centroids[i*2]));
		free(obj[0]);
		nb_container_destroy(area_trg);
		free(obj);
		i += 1;
	}
	return centroids;
}

static void _calculate_area_centroid(const vcn_mesh_t *mesh,
				    nb_container_t *area_trg,
				    double centroid[2])
{
	memset(centroid, 0, 2 * sizeof(*centroid));
	double* trg_centroids = 
		alloca(nb_container_get_length(area_trg) *
		       2 * sizeof(*trg_centroids));
	double vol = 0;
	uint32_t trg_count = 0;
	while (nb_container_is_not_empty(area_trg)) {
		msh_trg_t* trg = nb_container_delete_first(area_trg);
		/* Restore color */
		attr_t* trg_attr = trg->attr;
		trg->attr = trg_attr->data;
		free(trg_attr);
		/* Compute triangle centroid */
		trg_centroids[trg_count * 2] = (trg->v1->x[0] +
						trg->v2->x[0] +
						trg->v3->x[0]) /
			(3.0 * mesh->scale) + mesh->xdisp;
		trg_centroids[trg_count*2+1] = (trg->v1->x[1] +
						trg->v2->x[1] +
						trg->v3->x[1]) /
			(3.0 * mesh->scale) + mesh->ydisp;
		/* Compute triangle volume (area) */
		double trg_vol = vcn_utils2D_get_trg_area(trg->v1->x,
							  trg->v2->x,
							  trg->v3->x) /
			POW2(mesh->scale);
		/* Compute global centroid */
		centroid[0] += trg_vol * trg_centroids[trg_count * 2];
		centroid[1] += trg_vol * trg_centroids[trg_count*2+1];
		/* Compute global volume */
		vol += trg_vol;
		/* Increase triangle counter */
		trg_count += 1;
	}
	centroid[0] /= vol;
	centroid[1] /= vol;
	
	uint32_t idx_closest_centroid = 0;
	double min_dist = vcn_utils2D_get_dist2(trg_centroids, centroid);
	for (uint32_t i = 1; i < trg_count; i++) {
		double dist = 
			vcn_utils2D_get_dist2(&(trg_centroids[i*2]),
					      centroid);
		if (dist < min_dist) {
			min_dist = dist;
			idx_closest_centroid = i;
		}
	}
	memcpy(centroid, &(trg_centroids[idx_closest_centroid*2]),
	       2 * sizeof(*centroid));
}

static void calculate_area_centroid(const vcn_mesh_t *mesh,
				    nb_container_t *area_trg,
				    double centroid[2])
{
	msh_edge_t *max_edge = NULL;
	double max_length2 = 0;
	while (nb_container_is_not_empty(area_trg)) {
		msh_trg_t* trg = nb_container_delete_first(area_trg);
		/* Restore color */
		attr_t* trg_attr = trg->attr;
		trg->attr = trg_attr->data;
		free(trg_attr);
		
		max_edge = check_if_internal_edge_is_longer(trg->s1, max_edge,
							    &max_length2);
		max_edge = check_if_internal_edge_is_longer(trg->s2, max_edge,
							    &max_length2);
		max_edge = check_if_internal_edge_is_longer(trg->s3, max_edge,
							    &max_length2);
	}
	if (max_edge != NULL) {
		centroid[0] = (max_edge->v1->x[0] + max_edge->v2->x[0]) / 2.0;
		centroid[1] = (max_edge->v1->x[1] + max_edge->v2->x[1]) / 2.0;
		centroid[0] = centroid[0] / mesh->scale + mesh->xdisp;
		centroid[1] = centroid[1] / mesh->scale + mesh->ydisp;
	} else {
		memset(centroid, 0, 2 * sizeof(*centroid));
	}
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

static double colorize_infected_triangles(msh_trg_t* trg_infected,
					  nb_container_t* infected_trg,
					  bool blocking_with_input_segments)
{
	double area = 0.0;
	/* Colorize triangle */
	nb_container_insert(infected_trg, trg_infected);
	attr_t* trg_attr = calloc(1, sizeof(*trg_attr));
	trg_attr->id = 0x19FEC7ED;
	trg_attr->data = trg_infected->attr;
	trg_infected->attr = trg_attr;
	bool segment_nonblock = !blocking_with_input_segments;
	if (blocking_with_input_segments)
		segment_nonblock = !medge_is_subsgm(trg_infected->s1);
	if (trg_infected->t1 != NULL && segment_nonblock) {
		msh_trg_t* nb_trg = trg_infected->t1;
		if (NULL == nb_trg->attr) {
			area += colorize_infected_triangles(nb_trg,
							    infected_trg,
							    blocking_with_input_segments);
		} else {
			attr_t* nb_attr = (attr_t*)nb_trg->attr;
			if (nb_attr->id != 0x19FEC7ED)
				area += colorize_infected_triangles(nb_trg,
								    infected_trg,
								    blocking_with_input_segments);
		}
	}
	if (blocking_with_input_segments)
		segment_nonblock = !medge_is_subsgm(trg_infected->s2);
	if (NULL != trg_infected->t2 && segment_nonblock) {
		msh_trg_t* nb_trg = trg_infected->t2;
		if (NULL == nb_trg->attr) {
			area += colorize_infected_triangles(nb_trg,
							    infected_trg,
							    blocking_with_input_segments);
		} else {
			attr_t* nb_attr = nb_trg->attr;
			if (nb_attr->id != 0x19FEC7ED)
				area += colorize_infected_triangles(nb_trg,
								    infected_trg,
								    blocking_with_input_segments);
		}
	}
	if (blocking_with_input_segments)
		segment_nonblock = !medge_is_subsgm(trg_infected->s3);
	if (NULL != trg_infected->t3 && segment_nonblock) {
		msh_trg_t* nb_trg = trg_infected->t3;
		if (NULL == nb_trg->attr) {
			area += colorize_infected_triangles(nb_trg,
							    infected_trg,
							    blocking_with_input_segments);
		} else {
			attr_t* nb_attr = (attr_t*)nb_trg->attr;
			if (nb_attr->id != 0x19FEC7ED)
				area += colorize_infected_triangles(nb_trg,
								    infected_trg,
								    blocking_with_input_segments);
		}
	}

	area += vcn_utils2D_get_trg_area(trg_infected->v1->x,
					 trg_infected->v2->x,
					 trg_infected->v3->x);
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

double* vcn_mesh_get_centroids_of_enveloped_areas(const vcn_mesh_t *const mesh,
						  uint32_t* N_centroids)
{
	/* TEMPORAL (Repeated code): 
	 * Unify with code of vcn_mesh_get_centroid_of_subareas() */
	nb_container_t* areas = nb_container_create(NB_SORTED);
	nb_container_set_comparer(areas, compare_area1_isGreaterThan_area2);
	nb_iterator_t* iter = nb_iterator_create();
	nb_iterator_set_container(iter, mesh->ht_trg);
	while (nb_iterator_has_more(iter)) {
		msh_trg_t* trg = (msh_trg_t*) nb_iterator_get_next(iter);
		if (NULL != trg->attr) {
			attr_t* attr = trg->attr;
			if (attr->id == 0x19FEC7ED)
				continue;
		}
		nb_container_t* area_trg = nb_container_create(NB_SORTED);
		double* area = malloc(sizeof(*area));
		area[0] = colorize_infected_triangles(trg, area_trg, true);
		void** obj = malloc(2 * sizeof(*obj));
		obj[0] = area;
		obj[1] = area_trg;
		nb_container_insert(areas, obj);
	}
	nb_iterator_destroy(iter);
	/* Get centroids */
	N_centroids[0] = nb_container_get_length(areas);
	double *centroids = calloc(N_centroids[0] * 2, sizeof(*centroids));
	N_centroids[0] = 0;
	while (nb_container_is_not_empty(areas)) {
		void** obj = nb_container_delete_first(areas);
		free(obj[0]);
		nb_container_t* area_trg = obj[1];
		bool is_enveloped = true;
		double* trg_centroids = NULL;
		if (nb_container_is_not_empty(area_trg))
			trg_centroids = malloc(nb_container_get_length(area_trg) *
					       2 * sizeof(*trg_centroids));
		double* centroid = calloc(2, sizeof(*centroid));
		double vol = 0;
		uint32_t trg_count = 0;
		while (nb_container_is_not_empty(area_trg)) {
			msh_trg_t* trg = nb_container_delete_first(area_trg);
			/* Restore color */
			attr_t* trg_attr = (attr_t*)trg->attr;
			trg->attr = trg_attr->data;
			free(trg_attr);
			/* Check if it is enveloped */
			if (trg->t1 == NULL ||
			    trg->t2 == NULL || trg->t3 == NULL) {
				is_enveloped = false;
				break;
			}
			/* Compute triangle centroid */
			trg_centroids[trg_count * 2] = (trg->v1->x[0] +
							trg->v2->x[0] +
							trg->v3->x[0]) /
				(3.0 * mesh->scale) + mesh->xdisp;
			trg_centroids[trg_count*2+1] = (trg->v1->x[1] +
							trg->v2->x[1] +
							trg->v3->x[1]) /
				(3.0 * mesh->scale) + mesh->ydisp;
			/* Compute triangle volume (area) */
			double trg_vol = vcn_utils2D_get_trg_area(trg->v1->x,
								  trg->v2->x,
								  trg->v3->x) /
				POW2(mesh->scale);

			/* Compute global centroid */
			centroid[0] += trg_vol * trg_centroids[trg_count * 2];
			centroid[1] += trg_vol * trg_centroids[trg_count*2+1];

			/* Compute global volume */
			vol += trg_vol;

			/* Increase triangle counter */
			trg_count += 1;
		}

		/* Remove the attributes of the remaining triangles if
		 * it is enveloped */
		while (nb_container_is_not_empty(area_trg)) {
			msh_trg_t* trg = nb_container_delete_first(area_trg);
			/* Restore color */
			attr_t* trg_attr = (attr_t*)trg->attr;
			trg->attr = trg_attr->data;
			free(trg_attr);
		}

		if (is_enveloped) {
			centroid[0] /= vol;
			centroid[1] /= vol;
			uint32_t idx_closest_centroid = 0;
			double min_dist = vcn_utils2D_get_dist2(trg_centroids, centroid);
			for (uint32_t i = 1; i < trg_count; i++) {
				double dist = vcn_utils2D_get_dist2(&(trg_centroids[i*2]),
								    centroid);
				if (dist < min_dist) {
					min_dist = dist;
					idx_closest_centroid = i;
				}
			}
			centroids[N_centroids[0] * 2] =
				trg_centroids[idx_closest_centroid * 2];
			centroids[N_centroids[0]*2+1] =
				trg_centroids[idx_closest_centroid*2+1];
			N_centroids[0] += 1;
		}
		free(trg_centroids);
		free(centroid);
		nb_container_destroy(area_trg);
		free(obj);
	}
	/* Free memory */
	nb_container_destroy(areas);

	/* Output */
	double *aux = NULL;
	if (N_centroids[0] > 0) {
		aux = malloc(N_centroids[0] * 2 * sizeof(*aux));
		memcpy(aux, centroids, N_centroids[0] * 2 * sizeof(*aux));
	}
	free(centroids);
	return aux;
}

double vcn_mesh_clear_enveloped_areas(vcn_mesh_t* mesh,
				     double* area_removed)
{
	nb_container_t* areas = nb_container_create(NB_SORTED);
	nb_container_set_comparer(areas, compare_area1_isGreaterThan_area2);
	nb_iterator_t* iter = nb_iterator_create();
	nb_iterator_set_container(iter, mesh->ht_trg);
	while (nb_iterator_has_more(iter)) {
		msh_trg_t* trg = (msh_trg_t*) nb_iterator_get_next(iter);
		if (NULL != trg->attr) {
			attr_t* attr = trg->attr;
			if(attr->id == 0x19FEC7ED)
				continue;
		}
		nb_container_t* area_trg = nb_container_create(NB_SORTED);
		double* area = malloc(sizeof(*area));
		area[0] = colorize_infected_triangles(trg, area_trg, true);
		void** obj = malloc(2 * sizeof(*obj));
		obj[0] = area;
		obj[1] = area_trg;
		nb_container_insert(areas, obj);
	}
	nb_iterator_destroy(iter);

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
			if(NULL != ((attr_t*)trg->attr)->data)
				free(((attr_t*)trg->attr)->data);
			free(trg->attr);
			free(trg);
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
		attr_t* attr = (attr_t*)trg->attr;
		trg->attr = attr->data;
		free(attr);
	}
	nb_container_destroy(main_area);
	free(main_obj);
	/* Free memory */
	nb_container_destroy(areas);
  
	/* Rescale value */
	ret_val /= POW2(mesh->scale);
	if(NULL != area_removed)
		area_removed[0] /= POW2(mesh->scale);

	/* Return biggest area */
	return ret_val;
}

double vcn_mesh_keep_biggest_isolated_area(vcn_mesh_t* mesh,
				       double* area_removed)
{
	nb_container_t* areas = nb_container_create(NB_SORTED);
	nb_container_set_comparer(areas, compare_area1_isGreaterThan_area2);
	nb_iterator_t* iter = nb_iterator_create();
	nb_iterator_set_container(iter, mesh->ht_trg);
	while( nb_iterator_has_more(iter)) {
		msh_trg_t* trg = (msh_trg_t*)nb_iterator_get_next(iter);
		if (NULL != trg->attr) {
			attr_t* attr = trg->attr;
			if(attr->id == 0x19FEC7ED)
				continue;
		}
		nb_container_t* area_trg = nb_container_create(NB_SORTED);
		double* area = malloc(sizeof(*area));
		area[0] = colorize_infected_triangles(trg, area_trg, false);
		void** obj = malloc(2*sizeof(*obj));
		obj[0] = area;
		obj[1] = area_trg;
		nb_container_insert(areas, obj);
	}
	nb_iterator_destroy(iter);
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
			/* Restore color */
			attr_t* trg_attr = (attr_t*)trg->attr;
			if (NULL != trg_attr->data)
				free(trg_attr->data);
			free(trg_attr);
			free(trg);
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
		/* Restore color */
		attr_t* trg_attr = (attr_t*)trg->attr;
		trg->attr = trg_attr->data;
		free(trg_attr);
	}
	nb_container_destroy(main_area);
	free(main_obj);
	nb_container_destroy(areas);

	/* Rescale value */
	ret_val /= POW2(mesh->scale);
	if (NULL != area_removed)
		area_removed[0] /= POW2(mesh->scale);

	/* Return biggest area */
	return ret_val;
}

uint32_t vcn_mesh_delete_isolated_segments(vcn_mesh_t *const restrict mesh)
{
	uint32_t removed = 0;
	for (uint32_t i = 0; i < mesh->N_input_sgm; i++) {
		msh_edge_t* restrict sgm = mesh->input_sgm[i];
		if (NULL == sgm)
			continue;
		if (NULL != sgm->t1 || NULL != sgm->t2)
			continue;
		mesh->input_sgm[i] = NULL;
		while (sgm != NULL) {
			msh_edge_t* to_free = sgm;
			sgm = medge_subsgm_next(sgm);
			nb_container_delete(mesh->ht_edge, to_free);
			medge_destroy_subsgm_attribute(to_free);
			free(to_free);
			removed += 1;
		}
	}
	return removed;
}

uint32_t vcn_mesh_delete_internal_input_segments(vcn_mesh_t *const restrict mesh)
{
	uint32_t removed = 0;
	for (uint32_t i = 0; i < mesh->N_input_sgm; i++) {
		msh_edge_t* restrict sgm = mesh->input_sgm[i];
		if (NULL == sgm)
			continue;
		if (NULL == sgm->t1 || NULL == sgm->t2)
			continue;    
		mesh->input_sgm[i] = NULL;      
		while (NULL != sgm) {
			msh_edge_t* prev_sgm = sgm;
			sgm = medge_subsgm_next(sgm);
			medge_destroy_subsgm_attribute(prev_sgm);
			removed += 1;
		}
	}
	return removed;
}

uint32_t vcn_mesh_delete_isolated_vertices(vcn_mesh_t* mesh)
{
	nb_container_t* useful_vtx =
		alloca(nb_container_get_memsize(NB_QUEUE));
	nb_container_init(useful_vtx, NB_SORTED);

	nb_iterator_t* iter = alloca(nb_iterator_get_memsize());
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, mesh->ht_trg);
	while (nb_iterator_has_more(iter)) {
		msh_trg_t* trg = (msh_trg_t*)nb_iterator_get_next(iter);
		/* Track useful vertices */
		nb_container_insert(useful_vtx, trg->v1);
		nb_container_insert(useful_vtx, trg->v2);
		nb_container_insert(useful_vtx, trg->v3);
	}
	nb_iterator_finish(iter);

	nb_container_t *vtx_to_destroy = 
		alloca(nb_container_get_memsize(NB_QUEUE));
	nb_container_init(vtx_to_destroy, NB_QUEUE);
	
	uint32_t deleted = vcn_bins2D_get_length(mesh->ug_vtx);
	vcn_bins2D_iter_t* bins2D_iter = vcn_bins2D_iter_create();
	vcn_bins2D_iter_set_bins(bins2D_iter, mesh->ug_vtx);
	while (vcn_bins2D_iter_has_more(bins2D_iter)) {
		const msh_vtx_t* vtx = vcn_bins2D_iter_get_next(bins2D_iter);
		if(!nb_container_exist(useful_vtx, vtx)){
			vcn_bins2D_delete(mesh->ug_vtx, vtx);
			if (_NB_INPUT_VTX == vtx->attr ||
			    _NB_SUBSEGMENT_VTX == vtx->attr) {
				for (uint32_t i=0; i < mesh->N_input_vtx; i++) {
					if (vtx == mesh->input_vtx[i]) {
						mesh->input_vtx[i] = NULL;
						break;
					}
				}
			}
			nb_container_insert(vtx_to_destroy, vtx);
		}
	}
	vcn_bins2D_iter_destroy(bins2D_iter);
	nb_container_finish(useful_vtx);

	while (nb_container_is_not_empty(vtx_to_destroy)) {
		msh_vtx_t *vtx = nb_container_delete_first(vtx_to_destroy);
		free(vtx);
	}
	nb_container_finish(vtx_to_destroy);

	deleted = deleted - vcn_bins2D_get_length(mesh->ug_vtx);
	return deleted;
}
