#ifndef __NB_GEOMETRIC_BOT_KNN_BINS2D_H__
#define __NB_GEOMETRIC_BOT_KNN_BINS2D_H__

#include <stdint.h>
#include <stdbool.h>
#include "nb/container_bot/container.h"
#include "nb/geometric_bot/point2D.h"

typedef struct nb_bins2D_s nb_bins2D_t;

uint32_t nb_bins2D_get_memsize(void);
void nb_bins2D_init(nb_bins2D_t *bins2D, double size_of_bins);
void nb_bins2D_finish(nb_bins2D_t *bins2D);

nb_bins2D_t* nb_bins2D_create(double size_of_bins);
void nb_bins2D_destroy(nb_bins2D_t* bins2D);
void nb_bins2D_clear(nb_bins2D_t* bins2D);
void nb_bins2D_set_destroyer(nb_bins2D_t* bins2D,
			      void (*destroy)(void*));
void nb_bins2D_insert(nb_bins2D_t *const bins2D,
		       const nb_point2D_t *const point);
nb_point2D_t* nb_bins2D_delete(nb_bins2D_t* bins2D,
				 const nb_point2D_t *const point);
nb_point2D_t* nb_bins2D_delete_first(nb_bins2D_t* bins2D);

/**
 * @brief Calculate the m nearest neighbours available, with m <= k.
 * @return The output is sorted based on the dist(...) function.
 */
uint32_t nb_bins2D_get_knn(const nb_bins2D_t *const bins2D,
			    const nb_point2D_t *const p,
			    uint32_t k, nb_point2D_t* knn[],
			    double knn_dist[]);

void nb_bins2D_set_filter(nb_bins2D_t *bins2D, 
			   bool (*filter)(const nb_point2D_t *const p_ref,
					  const nb_point2D_t *const p,
					  const void *const data));

void nb_bins2D_set_filter_data(nb_bins2D_t *bins2D, const void *data);

void nb_bins2D_get_candidate_points_to_min_delaunay
				(const nb_bins2D_t *const bins2D,
				 const nb_point2D_t *const p1, 
				 const nb_point2D_t *const p2,
				 nb_container_t* vertices);

void nb_bins2D_get_points_inside_circle(const nb_bins2D_t *const bins2D,
					 double center[2],
					 double radius,
					 nb_container_t *points_inside);
bool nb_bins2D_are_points_inside_circle(const nb_bins2D_t *const bins2D,
					 double center[2],
					 double radius);
  
uint32_t nb_bins2D_get_N_bins(const nb_bins2D_t *const  bins2D);
uint32_t nb_bins2D_get_min_points_x_bin(const nb_bins2D_t *const bins2D);
uint32_t nb_bins2D_get_length(const nb_bins2D_t *const  bins2D);
bool nb_bins2D_is_empty(const nb_bins2D_t *const  bins2D);
bool nb_bins2D_is_not_empty(const nb_bins2D_t *const  bins2D);
double nb_bins2D_get_size_of_bins(const nb_bins2D_t *const bins2D);

#endif
