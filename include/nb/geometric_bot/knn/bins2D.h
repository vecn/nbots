#ifndef __NB_GEOMETRIC_BOT_KNN_BINS2D_H__
#define __NB_GEOMETRIC_BOT_KNN_BINS2D_H__

#include <stdint.h>
#include <stdbool.h>
#include "nb/container_bot/container.h"
#include "nb/geometric_bot/point2D.h"

typedef struct vcn_bins2D_s vcn_bins2D_t;

vcn_bins2D_t* vcn_bins2D_create(double size_of_bins);
void vcn_bins2D_destroy(vcn_bins2D_t* bins2D);
void vcn_bins2D_clear(vcn_bins2D_t* bins2D);
void vcn_bins2D_enable_point_destroyer(vcn_bins2D_t* bins2D);
void vcn_bins2D_disable_point_destroyer(vcn_bins2D_t* bins2D);
void vcn_bins2D_set_attribute_destroyer(vcn_bins2D_t* bins2D,
					void (*destroy)(void*));
void vcn_bins2D_insert(vcn_bins2D_t *const bins2D,
		       const vcn_point2D_t *const point);
vcn_point2D_t* vcn_bins2D_delete(vcn_bins2D_t* bins2D,
				 const vcn_point2D_t *const point);
vcn_point2D_t* vcn_bins2D_delete_first(vcn_bins2D_t* bins2D);

/**
 * @brief Calculate the m nearest neighbours available, with m <= k.
 * @return The output is sorted based on the dist(...) function.
 */
uint32_t vcn_bins2D_get_knn(const vcn_bins2D_t *const bins2D,
			    const vcn_point2D_t *const p,
			    uint32_t k, vcn_point2D_t* knn[],
			    double knn_dist[]);

void vcn_bins2D_set_filter(vcn_bins2D_t *bins2D, 
			   bool (*filter)(const vcn_point2D_t *const p_ref,
					  const vcn_point2D_t *const p,
					  const void *const data));

void vcn_bins2D_set_filter_data(vcn_bins2D_t *bins2D, const void *data);

nb_container_t* vcn_bins2D_get_candidate_points_to_min_delaunay
(const vcn_bins2D_t *const bins2D,
 const vcn_point2D_t *const p1, 
 const vcn_point2D_t *const p2);

nb_container_t* vcn_bins2D_get_points_inside_circle
(const vcn_bins2D_t *const bins2D,
 double center[2],
 double radius);
bool vcn_bins2D_are_points_inside_circle(const vcn_bins2D_t *const bins2D,
					 double center[2],
					 double radius);
  
uint32_t vcn_bins2D_get_N_bins(const vcn_bins2D_t *const  bins2D);
uint32_t vcn_bins2D_get_min_points_x_bin(const vcn_bins2D_t *const bins2D);
uint32_t vcn_bins2D_get_length(const vcn_bins2D_t *const  bins2D);
bool vcn_bins2D_is_empty(const vcn_bins2D_t *const  bins2D);
bool vcn_bins2D_is_not_empty(const vcn_bins2D_t *const  bins2D);
double vcn_bins2D_get_size_of_bins(const vcn_bins2D_t *const bins2D);

#endif
