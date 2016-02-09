#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "vcn/container_bot/container.h"
#include "vcn/container_bot/iterator.h"
#include "vcn/geometric_bot/point2D.h"
#include "vcn/geometric_bot/utils2D.h"
#include "vcn/geometric_bot/knn/bins2D.h"
#include "vcn/geometric_bot/knn/bins2D_iterator.h"

#include "bins2D_structs.h"

#define MAX(a,b) (((a)>(b))?(a):(b))
#define POW2(a) ((a)*(a))

typedef struct {
	double c[2];
	double r;
} circle_t;


static int get_bin_coord(double x, double bin_size);
static double get_bin_cartesian_coord(int x, double bin_size);
static uint32_t bin_hash_key(const void *const bin_ptr);
static bool bin_are_equal(const void* restrict bin1_ptr,
				 const void* restrict bin2_ptr);
static void bin_destroy(void *bin_ptr);
static void null_destroy(void *data);
static bool null_filter(const vcn_point2D_t *const p_ref,
			const vcn_point2D_t *const p,
			const void *const data);
static bin2D_t* get_bin(const vcn_container_t *const bins,
			int xcell, int ycell);

static uint32_t count_bins_in_block(const vcn_container_t *const bins,
				    const int block[4]);
static uint32_t count_bins_iterating_container
                              (const vcn_container_t *const restrict bins,
			       const int block[4]);
static uint32_t count_bins_iterating_block
                              (const vcn_container_t *const restrict bins,
			       const int block[4]);
static vcn_container_t* get_bins_block
                         (const vcn_container_t *const bins,
			  int add_block[4],
			  int substract_block[4]
			  /* NULL if not required */);
static bool bin_is_inside_block(int xcell, int ycell, const int block[4]);

static bool bin_is_outside_block(int xcell, int ycell, const int block[4]);

static vcn_container_t* knn_get_layer_bins
                     (const vcn_bins2D_t *const bins2D, int layer,
		      const vcn_point2D_t *const p);

static void knn_scan_list(const vcn_container_t *const points,
			  const vcn_point2D_t *const p_ref,
			  uint32_t k, vcn_point2D_t** knn,
			  uint32_t* m_nearest,
			  double* knn_dist,
			  bool (*filter)(const vcn_point2D_t *const p_ref,
					 const vcn_point2D_t *const p,
					 const void *const data),
			  const void* const filter_data);
static bool half_side_have_points(const vcn_bins2D_t* const bins2D,
				  const vcn_point2D_t *const p1,
				  const vcn_point2D_t *const p2);
static int8_t bin_which_half_side(const bin2D_t* const bin,
				  double bin_size,
				  const vcn_point2D_t *const p1,
				  const vcn_point2D_t *const p2);
static bool bin_have_points_in_half_side(const bin2D_t* const bin,
					 const vcn_point2D_t *const p1,
					 const vcn_point2D_t *const p2);
static vcn_container_t* delaunay_get_layer_bins
                     (const vcn_bins2D_t *const bins2D,
		      const vcn_point2D_t *const p1,
		      const vcn_point2D_t *const p2,
		      const circle_t *const layer_circle,
		      int layer);

static void set_block_from_circle(int block[4],
				  const vcn_bins2D_t *const bins2D,
				  const circle_t *const circle);
static void delaunay_get_points(const vcn_container_t *const points,
				const vcn_point2D_t *const p1,
				const vcn_point2D_t *const p2,
				const circle_t *const circle,
				vcn_container_t *vertices,
				vcn_container_t* outside_vtx);
static void delaunay_insert_if_is_candidate
				(vcn_container_t* vertices,
				vcn_container_t* outside_vtx,
				 const vcn_point2D_t *const restrict p1,
				 const vcn_point2D_t *const restrict p2,
				 const circle_t *const circle,
				 const vcn_point2D_t *const restrict p);
static bool is_in_half_side(const vcn_point2D_t *const p1,
			    const vcn_point2D_t *const p2,
			    const vcn_point2D_t *const p);
static bool is_inside_circle(const circle_t *const circle,
			     const vcn_point2D_t *const p);
static void set_circle_from_layer(circle_t *circle,
				  const vcn_point2D_t *const p1,
				  const vcn_point2D_t *const p2,
				  const vcn_bins2D_t *const bins2D,
				  int layer);
static void delaunay_set_inside_points_of_prev_layer
				(vcn_container_t* vertices,
				 vcn_container_t* outside_vtx,
				 const circle_t *const circle);

/****************** Functions implementation *************************/
static inline int get_bin_coord(double x, double bin_size)
{
	/* Shifting the values to avoid the double zero coordinates */
	return (int)(((x < 0)?(x - bin_size):x)/bin_size);
}

static inline double get_bin_cartesian_coord(int x, double bin_size)
{
	return bin_size * ((x < 0)?(x + 1):x);
}

static inline uint32_t bin_hash_key(const void *const bin_ptr)
{
	const bin2D_t *const bin = bin_ptr;
	return (uint32_t)
		((bin->x*73856093) ^ (bin->y*19349663) ^
		 (((bin->y != 0)?(bin->x/bin->y):(13))*83492791));
}

static inline bool bin_are_equal(const void* restrict bin1_ptr,
				 const void* restrict bin2_ptr)
{
	const bin2D_t *const restrict bin1 = bin1_ptr;
	const bin2D_t *const restrict bin2 = bin2_ptr;
	return (bin1->x == bin2->x) && (bin1->y == bin2->y);
}

static inline void bin_destroy(void *bin_ptr)
{
	bin2D_t* bin = bin_ptr;
	vcn_container_destroy(bin->points);
	free(bin);
}

vcn_bins2D_t* vcn_bins2D_create(double size_of_bins)
{
	vcn_bins2D_t* bins2D = calloc(1, sizeof(vcn_bins2D_t));
	bins2D->size_of_bins = size_of_bins;
	bins2D->bins = vcn_container_create(VCN_CONTAINER_HASH);
	vcn_container_set_key_generator(bins2D->bins, bin_hash_key);
	vcn_container_set_comparer(bins2D->bins, bin_are_equal);
	vcn_container_set_destroyer(bins2D->bins, bin_destroy);
	bins2D->destroy = (void(*)(void*)) vcn_point2D_destroy;
	bins2D->destroy_attr = null_destroy;
	bins2D->filter = null_filter;
	return bins2D;
}

static inline void null_destroy(void *data)
{
	; /* NULL statement */
}

static inline bool null_filter(const vcn_point2D_t *const p_ref,
			       const vcn_point2D_t *const p,
			       const void *const data)
{
	return true;
}

inline void vcn_bins2D_destroy(vcn_bins2D_t* bins2D)
{
	vcn_bins2D_clear(bins2D);
	vcn_container_destroy(bins2D->bins);
	free(bins2D);
}

void vcn_bins2D_clear(vcn_bins2D_t* bins2D)
{
	while (vcn_bins2D_is_not_empty(bins2D)) {
		vcn_point2D_t* point = vcn_bins2D_delete_first(bins2D);
		bins2D->destroy(point);
	}
	vcn_container_clear(bins2D->bins);
	bins2D->length = 0;
}

inline void vcn_bins2D_enable_point_destroyer(vcn_bins2D_t* bins2D)
{
	bins2D->destroy = (void (*)(void*)) vcn_point2D_destroy;
}

inline void vcn_bins2D_disable_point_destroyer(vcn_bins2D_t* bins2D)
{
	bins2D->destroy = null_destroy;
}

inline void vcn_bins2D_set_attribute_destroyer(vcn_bins2D_t* bins2D,
					       void (*destroy)(void*))
{
	bins2D->destroy_attr = destroy;
}

void vcn_bins2D_insert(vcn_bins2D_t *const restrict bins2D,
		       const vcn_point2D_t *const restrict point)
{
	bin2D_t key_bin;
	key_bin.x = get_bin_coord(point->x[0], bins2D->size_of_bins);
	key_bin.y = get_bin_coord(point->x[1], bins2D->size_of_bins);

	bin2D_t* bin = vcn_container_exist(bins2D->bins, &key_bin);
	if (NULL == bin) {
		bin = calloc(1, sizeof(bin2D_t));
		bin->x = key_bin.x;
		bin->y = key_bin.y;
		bin->points = vcn_container_create(VCN_CONTAINER_QUEUE);
		vcn_container_set_comparer(bin->points, vcn_point2D_are_equal);
		vcn_container_insert(bins2D->bins, bin);
	}
	vcn_container_insert(bin->points, point);
	bins2D->length += 1;
}

vcn_point2D_t* vcn_bins2D_delete(vcn_bins2D_t* bins2D,
				 const vcn_point2D_t *const point)
{
	bin2D_t key_bin;
	key_bin.x = get_bin_coord(point->x[0], bins2D->size_of_bins);
	key_bin.y = get_bin_coord(point->x[1], bins2D->size_of_bins);
	bin2D_t *bin =  vcn_container_exist(bins2D->bins, &key_bin);
	vcn_point2D_t *deleted_point = NULL;
	if (NULL != bin) {
		deleted_point = vcn_container_delete(bin->points, point);
		if (deleted_point == point) {
			bins2D->length -= 1;
			if (vcn_container_is_empty(bin->points)) {
				vcn_container_delete(bins2D->bins, bin);
				bin_destroy(bin);
			}
		}
	}
	return deleted_point;
}

vcn_point2D_t* vcn_bins2D_delete_first(vcn_bins2D_t* bins2D)
{
	vcn_point2D_t *point = NULL;
	bin2D_t *bin = vcn_container_get_first(bins2D->bins);
	if (NULL != bin) {
		point = vcn_container_delete_first(bin->points);
		if (vcn_container_is_empty(bin->points)) {
			vcn_container_delete(bins2D->bins, bin);
			bin_destroy(bin);
		}
	}
	return point;
}

static inline bin2D_t* get_bin(const vcn_container_t *const restrict bins,
			       int xcell, int ycell)
{
	bin2D_t key_cell;
	key_cell.x = xcell;
	key_cell.y = ycell;
	return vcn_container_exist(bins, &key_cell);
}

static inline uint32_t count_bins_in_block
                            (const vcn_container_t *const restrict bins,
			     const int block[4])
{
	uint32_t N_bins;
	if (vcn_container_get_length(bins) < 
	    (block[2] - block[0]) * (block[3] - block[1]))
		N_bins = count_bins_iterating_container(bins, block);
	else
		N_bins = count_bins_iterating_block(bins, block);
	return N_bins;
}

static uint32_t count_bins_iterating_container
                              (const vcn_container_t *const restrict bins,
			       const int block[4])
{
	uint32_t N_bins = 0;
	vcn_iterator_t* iter = vcn_iterator_create();
	vcn_iterator_set_container(iter, bins);
	while (vcn_iterator_has_more(iter)) {
		const bin2D_t* bin = vcn_iterator_get_next(iter);
		if (bin_is_inside_block(bin->x, bin->y, block))
			N_bins += 1;
	}
	vcn_iterator_destroy(iter);
	return N_bins;
}

static uint32_t count_bins_iterating_block
                              (const vcn_container_t *const restrict bins,
			       const int block[4])
{
	uint32_t N_bins = 0;
	for (int i = block[0]; i <= block[2]; i++) {
		for (int j = block[1]; j <= block[3]; j++) {
			if (NULL != get_bin(bins, i, j))
				N_bins += 1;
		}
	}
	return N_bins;
}

static inline vcn_container_t* get_bins_block
                             (const vcn_container_t *const restrict bins,
			      int add_block[4],
			      int substract_block[4]
			      /* NULL if not required */)
{
	vcn_container_t* bins_block = vcn_container_create(VCN_CONTAINER_QUEUE);
	if (vcn_container_get_length(bins) < 
	    (add_block[2] - add_block[0]) * (add_block[3] - add_block[1])) {
		/* Iterate cells in the hash table */
		vcn_iterator_t* iter = vcn_iterator_create();
		vcn_iterator_set_container(iter, bins);
		while (vcn_iterator_has_more(iter)) {
			const bin2D_t* bin = vcn_iterator_get_next(iter);
			if (bin_is_outside_block(bin->x, bin->y, add_block))
				continue;
			if (NULL != substract_block)
				if (bin_is_inside_block(bin->x, bin->y, substract_block)) 
					continue;
			/* Insert points */
			vcn_container_insert(bins_block, bin->points);	
		}
		vcn_iterator_destroy(iter);
	} else {
		/* Iterate cells in the block */
		for (int i = add_block[0]; i <= add_block[2]; i++) {
			for (int j = add_block[1]; j <= add_block[3]; j++) {
				if (NULL != substract_block)
					if (bin_is_inside_block(i, j, substract_block)) 
						continue;
				bin2D_t *bin = get_bin(bins, i, j);
				if (NULL == bin)
					continue;
				/* Insert points */
				vcn_container_insert(bins_block, bin->points);
			}
		}
	}
	return bins_block;
}

static inline bool bin_is_inside_block(int xcell, int ycell,
				       const int block[4])
{
	bool inside_xspan = (xcell >= block[0] && xcell <= block[2]);
	bool inside_yspan = (ycell >= block[1] && ycell <= block[3]);
	return inside_xspan && inside_yspan;
}

static inline bool bin_is_outside_block(int xcell, int ycell,
					const int block[4])
{
	bool outside_xspan = (xcell < block[0] || xcell > block[2]);
	bool outside_yspan = (ycell < block[1] || ycell > block[3]);
	return outside_xspan || outside_yspan;
}

static inline vcn_container_t* knn_get_layer_bins
                     (const vcn_bins2D_t *const restrict bins2D, 
		      int layer,
		      const vcn_point2D_t *const restrict p)
{
	/* The zero-layer is conformed  by the bin containing the point of 
	 * reference.
	 * To calculate the subsequent layers,  we must increase the radius
	 * of an initial circle until obtain a bigger enveloping rectangle.
	 * A Layer does not conatin any bin of the previous layers.
	 */
  
	/* Compute spatial-bin coordinates of the point */
	int xbin = get_bin_coord(p->x[0], bins2D->size_of_bins);
	int ybin = get_bin_coord(p->x[1], bins2D->size_of_bins);

	/* Create list to store the bins of the layer */
	vcn_container_t* layer_bins = vcn_container_create(VCN_CONTAINER_QUEUE);

	if (0 == layer) {
		bin2D_t* bin = get_bin(bins2D->bins, xbin, ybin);
		if (NULL != bin)
			/* Insert points */
			vcn_container_insert(layer_bins, bin->points);
	} else if(0 < layer) {
		/* Get bottom horizontal bins */
		for (int i = xbin - layer; i <= xbin + layer; i++) {
			bin2D_t* bin = get_bin(bins2D->bins, i, ybin - layer);
			if (NULL != bin)
				vcn_container_insert(layer_bins, bin->points);
		}
		/* Get top horizontal bins */
		for (int i = xbin - layer; i <= xbin + layer; i++) {
			bin2D_t* bin = get_bin(bins2D->bins, i, ybin + layer);
			if (NULL != bin)
				vcn_container_insert(layer_bins, bin->points);
		}
		/* Get left vertical bins */
		for (int j = ybin - layer + 1; j < ybin + layer; j++) {
			bin2D_t* bin = get_bin(bins2D->bins, xbin - layer, j);
			if (NULL != bin)
				vcn_container_insert(layer_bins, bin->points);
		}
		/* Get right vertical bins */
		for (int j = ybin - layer + 1; j < ybin + layer; j++) {
			bin2D_t* bin = get_bin(bins2D->bins, xbin + layer, j);
			if (NULL != bin)
				vcn_container_insert(layer_bins, bin->points);
		}
	}
	return layer_bins;
}

static inline void knn_scan_list(const vcn_container_t *const points,
				 const vcn_point2D_t *const restrict p_ref,
				 uint32_t k, vcn_point2D_t** knn,
				 uint32_t* m_nearest,
				 double* knn_dist,
				 bool (*filter)(const vcn_point2D_t *const p_ref,
						const vcn_point2D_t *const p,
						const void *const data),
				 const void* const restrict filter_data)
{
	/* Iterate points in the cell */
	vcn_iterator_t *const restrict iter = vcn_iterator_create();
	vcn_iterator_set_container(iter, points);
	while (vcn_iterator_has_more(iter)) {
		/* Get next point */
		vcn_point2D_t *restrict p =
			(vcn_point2D_t*) vcn_iterator_get_next(iter);

		/* Filter if the pointer is the same */
		if (p_ref == p)
			continue;
    
		/* Filter point */
		if(!filter(p_ref, p, filter_data)) 
			continue;
    
		/* Calculate distance */
		double dist_p = vcn_utils2D_get_dist(p_ref->x, p->x);

		for (uint32_t i = 0; i < m_nearest[0]; i++) {
			/* Minimize distance inserting in ascending order */
			if (dist_p < knn_dist[i]) {
				/* Swap value */
				double aux = dist_p;
				dist_p = knn_dist[i];
				knn_dist[i] = aux;
				vcn_point2D_t* p_aux = p;
				p = knn[i];
				knn[i] = p_aux;
			}
		}
		/* Add new nearest neighbour to the output if there are space */
		if (m_nearest[0] < k) {
			knn_dist[m_nearest[0]] = dist_p;
			knn[m_nearest[0]] = p;
			m_nearest[0] += 1;
		}
	}
	vcn_iterator_destroy(iter);
}

uint32_t vcn_bins2D_get_knn(const vcn_bins2D_t *const restrict bins2D,
			    const vcn_point2D_t *const restrict p,
			    uint32_t k, vcn_point2D_t* knn[], double knn_dist[])
/* Returns the m nearest neighbours available, with m <= k.
 * The output is sorted based on the dist(...) function.
 */
{  
	/* Scan bins layer by layer. In the following example, the dot
	 * represents the centroid of the reference, and the labels are
	 * the number of the layer containing the cell.
	 * This is the layer distribution using euclidean distance:
	 *
	 *                                |2|2|2|2|2| <--- Bins
	 *                                |2|1|1|1|2|
	 *       Spatial bin Layers --->  |2|1|·|1|2|
	 *                                |2|1|1|1|2|
	 *                                |2|2|2|2|2|
	 */
	uint32_t n_scanned = 0;
	uint32_t m_nearest = 0;
	int layer = 0;
	while (n_scanned < bins2D->length  && m_nearest < k) {
		/* 1. Get cells where the search should be performed. */
		vcn_container_t* layer_cells = knn_get_layer_bins(bins2D, layer, p);

		/* 2. Compute nearest neighbours */
		while (vcn_container_is_not_empty(layer_cells)) {
			const vcn_container_t *const restrict points =
				vcn_container_delete_first(layer_cells);
			n_scanned += vcn_container_get_length(points);
			knn_scan_list(points, p, k, knn, &m_nearest, knn_dist, 
				      bins2D->filter,
				      bins2D->filter_data);
		}
		/* Destroy list that had the layer cells */
		vcn_container_destroy(layer_cells);

		/* Increment layer */
		layer += 1;
	}

	/* Return number of closest points */
	return m_nearest;
}

inline void vcn_bins2D_set_filter(vcn_bins2D_t *bins2D, 
				  bool (*filter)
				  	(const vcn_point2D_t *const p_ref,
					 const vcn_point2D_t *const p,
					 const void *const data))
{
	bins2D->filter = filter;
}

inline void vcn_bins2D_set_filter_data(vcn_bins2D_t *bins2D, const void *data)
{
	bins2D->filter_data = data;
}

static bool half_side_have_points(const vcn_bins2D_t* const bins2D,
				  const vcn_point2D_t *const p1,
				  const vcn_point2D_t *const p2)
{	
	vcn_iterator_t* iter = vcn_iterator_create();
	vcn_iterator_set_container(iter, bins2D->bins);
	bool have_points = false;
	while (vcn_iterator_has_more(iter)) {
		const bin2D_t* bin = vcn_iterator_get_next(iter);
		int8_t half_side = 
			bin_which_half_side(bin, bins2D->size_of_bins, p1, p2);
		if (1 == half_side && vcn_container_is_not_empty(bin->points)) {
			have_points = true;
			break;
		} else if (0 == half_side) {
			if (bin_have_points_in_half_side(bin, p1, p2)) {
				have_points = true;
				break;
			}
		}
	}
	vcn_iterator_destroy(iter);
	return have_points;
}

static inline int8_t bin_which_half_side(const bin2D_t* const bin,
					 double bin_size,
					 const vcn_point2D_t *const p1,
					 const vcn_point2D_t *const p2)
{
	vcn_point2D_t bin_p;
	bin_p.x[0] = get_bin_cartesian_coord(bin->x, bin_size);
	bin_p.x[1] = get_bin_cartesian_coord(bin->y, bin_size);
	double cp[2];
	vcn_utils2D_get_closest_pnt_to_sgm(p1->x, p2->x,
					   bin_p.x, cp);
	double dist2 = vcn_utils2D_get_dist2(bin_p.x, cp);
	int8_t half_side;
	if (dist2 < POW2(bin_size)) {
		half_side = 0;
	} else {
		if (is_in_half_side(p1, p2, &bin_p))
			half_side = 1;
		else
			half_side = -1;
	}
	return half_side;
}

static bool bin_have_points_in_half_side(const bin2D_t* const bin,
					 const vcn_point2D_t *const p1,
					 const vcn_point2D_t *const p2)
{
	vcn_iterator_t *iter = vcn_iterator_create();
	vcn_iterator_set_container(iter, bin->points);
	bool have_points = false;
	while (vcn_iterator_has_more(iter)) {
		const vcn_point2D_t *p = vcn_iterator_get_next(iter);
		if (is_in_half_side(p1, p2, p)) {
			have_points = true;
			break;
		}
	}
	vcn_iterator_destroy(iter);
	return have_points;
}

static vcn_container_t* delaunay_get_layer_bins
			(const vcn_bins2D_t *const bins2D,
			 const vcn_point2D_t *const p1,
			 const vcn_point2D_t *const p2,
			 const circle_t *const layer_circle,
			 int layer)
{
	int layer_block[4];
	set_block_from_circle(layer_block, bins2D, layer_circle);

	int *substract_block = NULL;
	int interior_block[4];
	if (0 < layer) {
		circle_t c_in;
		set_circle_from_layer(&c_in, p1, p2, bins2D, layer - 1);
		set_block_from_circle(interior_block, bins2D, &c_in);
		substract_block = interior_block;
	}
	return get_bins_block(bins2D->bins, layer_block, substract_block);
}

static void set_block_from_circle(int block[4],
				  const vcn_bins2D_t *const bins2D,
				  const circle_t *const circle)
{
	block[0] = get_bin_coord(circle->c[0] - circle->r, bins2D->size_of_bins);
	block[1] = get_bin_coord(circle->c[1] - circle->r, bins2D->size_of_bins);
	block[2] = get_bin_coord(circle->c[0] + circle->r, bins2D->size_of_bins);
	block[3] = get_bin_coord(circle->c[1] + circle->r, bins2D->size_of_bins);

}

static void delaunay_get_points(const vcn_container_t *const restrict points,
				const vcn_point2D_t *const restrict p1,
				const vcn_point2D_t *const restrict p2,
				const circle_t *const circle,
				vcn_container_t* vertices,
				vcn_container_t* outside_vtx)
{
	vcn_iterator_t *const restrict iter = vcn_iterator_create();
	vcn_iterator_set_container(iter, points);
	while (vcn_iterator_has_more(iter)) {
		const vcn_point2D_t* restrict p = vcn_iterator_get_next(iter);
		delaunay_insert_if_is_candidate(vertices, outside_vtx,
						p1, p2, circle, p);
	}
	vcn_iterator_destroy(iter);
}

static void delaunay_insert_if_is_candidate
				(vcn_container_t* vertices,
				 vcn_container_t* outside_vtx,
				 const vcn_point2D_t *const restrict p1,
				 const vcn_point2D_t *const restrict p2,
				 const circle_t *const circle,
				 const vcn_point2D_t *const restrict p)
{
	if (p1 != p && p2 != p) {
		if (is_in_half_side(p1, p2, p)) {
			if (is_inside_circle(circle, p))
				vcn_container_insert(vertices, p);
			else
				vcn_container_insert(outside_vtx, p);
		}
	}
}

static inline bool is_in_half_side(const vcn_point2D_t *const restrict p1,
				   const vcn_point2D_t *const restrict p2,
				   const vcn_point2D_t *const restrict p)
{
	double sign = vcn_utils2D_get_2x_trg_area(p1->x, p2->x, p->x);
	return (sign >= VCN_GEOMETRIC_TOL);
}

static inline bool is_inside_circle(const circle_t *const circle,
				    const vcn_point2D_t *const p)
{
	double dist2 = vcn_utils2D_get_dist2(circle->c, p->x);
	return (POW2(circle->r) > dist2);
}

vcn_container_t* vcn_bins2D_get_candidate_points_to_min_delaunay
                   (const vcn_bins2D_t *const restrict bins2D,
		    const vcn_point2D_t *const restrict p1, 
		    const vcn_point2D_t *const restrict p2)
{
	/* Scan bins layer by layer. In the following example, the dots
	 * represents the vertices of the segment, and the labels are
	 * the number of the layer containing the bin.
	 * This is the layer distribution using delaunay distance:
	 *                              
	 *                                |3|3|3|3|3| <--- Layers
	 *                                |2|2|2|2|3|
	 *                                |·|1|1|2|3|
	 *                     Segment ---|>|\|1|2|3|
	 *                                | | |·|2|3|
	 */
	vcn_container_t *vertices = vcn_container_create(VCN_CONTAINER_QUEUE);
	vcn_container_t *outside_vtx = vcn_container_create(VCN_CONTAINER_QUEUE);
	bool have_points = half_side_have_points(bins2D, p1, p2);
	int layer = 0;
	while (vcn_container_is_empty(vertices) && have_points) {
		circle_t circle;
		set_circle_from_layer(&circle, p1, p2, bins2D, layer);

		delaunay_set_inside_points_of_prev_layer(vertices, outside_vtx,
							 &circle);

		/* 1. Get bins where the search should be performed. */
		vcn_container_t* layer_bins = 
			delaunay_get_layer_bins(bins2D, p1, p2, &circle, layer);

		/* 2. Get points from the bins of the layer */
		while (vcn_container_is_not_empty(layer_bins)) {
			const vcn_container_t *const restrict points =
				vcn_container_delete_first(layer_bins);
			delaunay_get_points(points, p1, p2, &circle,
					    vertices, outside_vtx);
		}
		vcn_container_destroy(layer_bins);	
		layer += 1;
	}
	vcn_container_destroy(outside_vtx);
	return vertices;
}

static void set_circle_from_layer(circle_t *circle,
				  const vcn_point2D_t *const p1,
				  const vcn_point2D_t *const p2,
				  const vcn_bins2D_t *const bins2D,
				  int layer)
{
	double n[2];
	n[0] = p1->x[1] - p2->x[1];
	n[1] = p2->x[0] - p1->x[0];

	double dist = sqrt(POW2(n[0]) + POW2(n[1]));
	n[0] /= dist;
	n[1] /= dist;

	circle->r = dist + layer * MAX(dist, bins2D->size_of_bins);
	double a = sqrt(POW2(circle->r) - 0.25 * POW2(dist));
	circle->c[0] = 0.5 * (p1->x[0] + p2->x[0]) + a * n[0];
	circle->c[1] = 0.5 * (p1->x[1] + p2->x[1]) + a * n[1];
}

static void delaunay_set_inside_points_of_prev_layer
				(vcn_container_t* vertices,
				 vcn_container_t* outside_vtx,
				 const circle_t *const circle)
{
	int container_id = vcn_container_get_id(outside_vtx);
	vcn_container_t *aux = vcn_container_create(container_id);
	while (vcn_container_is_not_empty(outside_vtx)) {
		vcn_point2D_t *p = vcn_container_delete_first(outside_vtx);
		if (is_inside_circle(circle, p))
			vcn_container_insert(vertices, p);
		else
			vcn_container_insert(aux, p);
	}
	vcn_container_merge(outside_vtx, aux);
	vcn_container_destroy(aux);
}

vcn_container_t* vcn_bins2D_get_points_inside_circle
					(const vcn_bins2D_t *const bins2D,
					 double center[2],
					 double radius)
{
	vcn_container_t* output = vcn_container_create(VCN_CONTAINER_QUEUE);
	/* Get block */
	int block[4];
	block[0] = get_bin_coord(center[0] - radius, bins2D->size_of_bins);
	block[1] = get_bin_coord(center[1] - radius, bins2D->size_of_bins);
	block[2] = get_bin_coord(center[0] + radius, bins2D->size_of_bins);
	block[3] = get_bin_coord(center[1] + radius, bins2D->size_of_bins);

	vcn_container_t* bins_block = get_bins_block(bins2D->bins, block, NULL);

	/* Compute points inside the circle */
	while (vcn_container_is_not_empty(bins_block)) {
		vcn_container_t* points = vcn_container_delete_first(bins_block);
		vcn_iterator_t* iter = vcn_iterator_create();
		vcn_iterator_set_container(iter, points);
		while (vcn_iterator_has_more(iter)) {
			const vcn_point2D_t *p = vcn_iterator_get_next(iter);
			double dist2 = vcn_utils2D_get_dist2(p->x, center);
			if (POW2(radius) - dist2 > VCN_GEOMETRIC_TOL_POW2)
				vcn_container_insert(output, p);
		}
		vcn_iterator_destroy(iter);
	}
	vcn_container_destroy(bins_block);

	return output;
}

bool vcn_bins2D_are_points_inside_circle
                        (const vcn_bins2D_t *const bins2D,
			 double center[2], double radius)
{
	/* Get block */
	int block[4];
	block[0]  = get_bin_coord(center[0] - radius, bins2D->size_of_bins);
	block[1]  = get_bin_coord(center[1] - radius, bins2D->size_of_bins);
	block[2]  = get_bin_coord(center[0] + radius, bins2D->size_of_bins);
	block[3]  = get_bin_coord(center[1] + radius, bins2D->size_of_bins);

	for (int i = block[0]; i <= block[2]; i++) {
		for (int j = block[1]; j <= block[3]; j++) {
			bin2D_t *bin = get_bin(bins2D->bins, i, j);
			if (NULL == bin)
				continue;
			vcn_container_t* points = bin->points;
			vcn_iterator_t* iter = vcn_iterator_create();
			vcn_iterator_set_container(iter, points);
			while (vcn_iterator_has_more(iter)) {
				const vcn_point2D_t *p = vcn_iterator_get_next(iter);
				double dist2 = vcn_utils2D_get_dist2(p->x, center);
				if (POW2(radius) - dist2 > VCN_GEOMETRIC_TOL_POW2) {
					vcn_iterator_destroy(iter);
					return true;
				}
			}
			vcn_iterator_destroy(iter);
		}
	}
	return false;
}

inline uint32_t vcn_bins2D_get_N_bins(const vcn_bins2D_t *const restrict bins2D)
{
	return vcn_container_get_length(bins2D->bins);
}

inline bool vcn_bins2D_is_empty(const vcn_bins2D_t *const restrict bins2D)
{
	return vcn_container_is_empty(bins2D->bins);
}

inline bool vcn_bins2D_is_not_empty(const vcn_bins2D_t *const restrict bins2D)
{
	return vcn_container_is_not_empty(bins2D->bins);
}

uint32_t vcn_bins2D_get_min_points_x_bin(const vcn_bins2D_t *const bins2D)
{
	uint32_t min = bins2D->length;
	vcn_iterator_t* iter = vcn_iterator_create();
	vcn_iterator_set_container(iter, bins2D->bins);
	while (vcn_iterator_has_more(iter)) {
		const bin2D_t *bin = vcn_iterator_get_next(iter);
		if (min > vcn_container_get_length(bin->points))
			min = vcn_container_get_length(bin->points);
	}
	vcn_iterator_destroy(iter);
	return min;
}

inline uint32_t vcn_bins2D_get_length(const vcn_bins2D_t *const restrict bins2D)
{
	return bins2D->length;
}

inline double vcn_bins2D_get_size_of_bins
				(const vcn_bins2D_t *const restrict bins2D)
{
	return bins2D->size_of_bins;
}