#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "nb/memory_bot.h"
#include "nb/container_bot/container.h"
#include "nb/container_bot/iterator.h"
#include "nb/geometric_bot/point2D.h"
#include "nb/geometric_bot/utils2D.h"
#include "nb/geometric_bot/knn/bins2D.h"

#include "bins2D_structs.h"

#define MAX(a,b) (((a)>(b))?(a):(b))
#define POW2(a) ((a)*(a))

typedef struct {
	double c[2];
	double r;
} circle_t;

static int get_bin_coord(double x, double bin_size);
static uint32_t bin_hash_key(const void *const bin_ptr);
static int8_t bin_compare(const void* restrict bin1_ptr,
			  const void* restrict bin2_ptr);
static void bin_destroy(void *bin_ptr);
static bool null_filter(const nb_point2D_t *const p_ref,
			const nb_point2D_t *const p,
			const void *const data);
static void destroy_points(nb_bins2D_t* bins2D);
static bin2D_t* get_bin(const nb_container_t *const bins,
			int xcell, int ycell);

static uint32_t count_bins_in_block(const nb_container_t *const bins,
				    const int block[4]);
static uint32_t count_bins_iterating_container
                              (const nb_container_t *const restrict bins,
			       const int block[4]);
static uint32_t count_bins_iterating_block
                              (const nb_container_t *const restrict bins,
			       const int block[4]);
static void get_bins_block(nb_container_t* bins_block,
			   const nb_container_t *const bins,
			   int add_block[4],
			   int substract_block[4]
			   /* NULL if not required */);
static void iterate_cells_in_container(nb_container_t* bins_block,
				       const nb_container_t *bins,
				       int add_block[4],
				       int substract_block[4]
				       /* NULL if not required */);
static void iterate_cells_in_block(nb_container_t* bins_block,
				   const nb_container_t *bins,
				   int add_block[4],
				   int substract_block[4]
				   /* NULL if not required */);
static bool bin_is_inside_block(int xcell, int ycell, const int block[4]);

static bool bin_is_outside_block(int xcell, int ycell, const int block[4]);

static void knn_get_layer_bins(const nb_bins2D_t *const bins2D, int layer,
			       const nb_point2D_t *const p,
			       nb_container_t *layer_bins);

static void knn_scan_list(const nb_container_t *const points,
			  const nb_point2D_t *const p_ref,
			  uint32_t k, nb_point2D_t** knn,
			  uint32_t* m_nearest,
			  double* knn_dist,
			  bool (*filter)(const nb_point2D_t *const p_ref,
					 const nb_point2D_t *const p,
					 const void *const data),
			  const void* const filter_data);
static bool half_side_have_points(const nb_bins2D_t* const bins2D,
				  const nb_point2D_t *const p1,
				  const nb_point2D_t *const p2);
static int8_t bin_which_half_side(const bin2D_t* const bin,
				  double bin_size,
				  const nb_point2D_t *const p1,
				  const nb_point2D_t *const p2);
static bool bin_have_points_in_half_side(const bin2D_t* const bin,
					 const nb_point2D_t *const p1,
					 const nb_point2D_t *const p2);
static void delaunay_get_layer_bins
                     (nb_container_t* layer_bins,
		      const nb_bins2D_t *const bins2D,
		      const nb_point2D_t *const p1,
		      const nb_point2D_t *const p2,
		      const circle_t *const layer_circle,
		      int layer);

static void set_block_from_circle(int block[4],
				  const nb_bins2D_t *const bins2D,
				  const circle_t *const circle);
static void delaunay_get_points(const nb_container_t *const points,
				const nb_point2D_t *const p1,
				const nb_point2D_t *const p2,
				const circle_t *const circle,
				nb_container_t *vertices,
				nb_container_t* outside_vtx);
static void delaunay_insert_if_is_candidate
				(nb_container_t* vertices,
				nb_container_t* outside_vtx,
				 const nb_point2D_t *const restrict p1,
				 const nb_point2D_t *const restrict p2,
				 const circle_t *const circle,
				 const nb_point2D_t *const restrict p);
static bool is_inside_circle(const circle_t *const circle,
			     const nb_point2D_t *const p);
static void set_circle_from_layer(circle_t *circle,
				  const nb_point2D_t *const p1,
				  const nb_point2D_t *const p2,
				  const nb_bins2D_t *const bins2D,
				  int layer);
static void delaunay_set_inside_points_of_prev_layer
				(nb_container_t* vertices,
				 nb_container_t* outside_vtx,
				 const circle_t *const circle);
static void get_points_inside_circle(nb_container_t *bins_block,
				     double center[2], double radius,
				     nb_container_t *points_inside);
static bool are_points_inside_circle(const nb_container_t *points,
				     double center[2], double radius);
/****************** Functions implementation *************************/
static inline int get_bin_coord(double x, double bin_size)
{
	/* Shifting the values to avoid the double zero coordinates */
	return (int)(((x < 0)?(x - bin_size):x)/bin_size);
}

static inline uint32_t bin_hash_key(const void *const bin_ptr)
{
	const bin2D_t *const bin = bin_ptr;
	return (uint32_t)
		((bin->x*73856093) ^ (bin->y*19349663) ^
		 (((bin->y != 0)?(bin->x/bin->y):(13))*83492791));
}

static inline int8_t bin_compare(const void* restrict bin1_ptr,
				 const void* restrict bin2_ptr)
{
	const bin2D_t *const restrict bin1 = bin1_ptr;
	const bin2D_t *const restrict bin2 = bin2_ptr;
	return (bin1->x == bin2->x) && (bin1->y == bin2->y) ? 0:1;
}

static inline void bin_destroy(void *bin_ptr)
{
	bin2D_t* bin = bin_ptr;
	nb_container_destroy(bin->points);
	nb_free_mem(bin);
}

uint32_t nb_bins2D_get_memsize(void)
{
	return sizeof(nb_bins2D_t);
}

void nb_bins2D_init(nb_bins2D_t *bins2D, double size_of_bins)
{
	memset(bins2D, 0, nb_bins2D_get_memsize());
	bins2D->size_of_bins = size_of_bins;
	bins2D->bins = nb_container_create(NB_HASH);
	nb_container_set_key_generator(bins2D->bins, bin_hash_key);
	nb_container_set_comparer(bins2D->bins, bin_compare);
	nb_container_set_destroyer(bins2D->bins, bin_destroy);
	bins2D->destroy = NULL;
	bins2D->filter = null_filter;
}

void nb_bins2D_finish(nb_bins2D_t *bins2D)
{
	destroy_points(bins2D);
	nb_container_destroy(bins2D->bins);
}

nb_bins2D_t* nb_bins2D_create(double size_of_bins)
{
	nb_bins2D_t* bins2D = nb_allocate_mem(nb_bins2D_get_memsize());
	nb_bins2D_init(bins2D, size_of_bins);
	return bins2D;
}

static inline bool null_filter(const nb_point2D_t *const p_ref,
			       const nb_point2D_t *const p,
			       const void *const data)
{
	return true;
}

void nb_bins2D_destroy(nb_bins2D_t* bins2D)
{
	nb_bins2D_finish(bins2D);
	nb_free_mem(bins2D);
}

void nb_bins2D_clear(nb_bins2D_t* bins2D)
{
	destroy_points(bins2D);
	nb_container_clear(bins2D->bins);
	bins2D->length = 0;
}

static void destroy_points(nb_bins2D_t* bins2D)
{
	if (NULL != bins2D->destroy) {
		while (nb_bins2D_is_not_empty(bins2D)) {
			nb_point2D_t* point = nb_bins2D_delete_first(bins2D);
			bins2D->destroy(point);
		}
	}
}

void nb_bins2D_set_destroyer(nb_bins2D_t* bins2D,
			      void (*destroy)(void*))
{
	bins2D->destroy = destroy;

}

void nb_bins2D_insert(nb_bins2D_t *const restrict bins2D,
		       const nb_point2D_t *const restrict point)
{
	bin2D_t key_bin;
	key_bin.x = get_bin_coord(point->x[0], bins2D->size_of_bins);
	key_bin.y = get_bin_coord(point->x[1], bins2D->size_of_bins);

	bin2D_t* bin = nb_container_exist(bins2D->bins, &key_bin);
	if (NULL == bin) {
		bin = nb_allocate_zero_mem(sizeof(bin2D_t));
		bin->x = key_bin.x;
		bin->y = key_bin.y;
		bin->points = nb_container_create(NB_QUEUE);
		nb_container_set_comparer(bin->points, nb_point2D_compare);
		nb_container_insert(bins2D->bins, bin);
	}
	nb_container_insert(bin->points, point);
	bins2D->length += 1;
}

nb_point2D_t* nb_bins2D_delete(nb_bins2D_t* bins2D,
				 const nb_point2D_t *const point)
{
	bin2D_t key_bin;
	key_bin.x = get_bin_coord(point->x[0], bins2D->size_of_bins);
	key_bin.y = get_bin_coord(point->x[1], bins2D->size_of_bins);
	bin2D_t *bin =  nb_container_exist(bins2D->bins, &key_bin);
	nb_point2D_t *deleted_point = NULL;
	if (NULL != bin) {
		deleted_point = nb_container_delete(bin->points, point);
		if (deleted_point == point) {
			bins2D->length -= 1;
			if (nb_container_is_empty(bin->points)) {
				nb_container_delete(bins2D->bins, bin);
				bin_destroy(bin);
			}
		}
	}
	return deleted_point;
}

nb_point2D_t* nb_bins2D_delete_first(nb_bins2D_t* bins2D)
{
	nb_point2D_t *point = NULL;
	bin2D_t *bin = nb_container_get_first(bins2D->bins);
	if (NULL != bin) {
		point = nb_container_delete_first(bin->points);
		if (nb_container_is_empty(bin->points)) {
			nb_container_delete(bins2D->bins, bin);
			bin_destroy(bin);
		}
	}
	return point;
}

static inline bin2D_t* get_bin(const nb_container_t *const restrict bins,
			       int xcell, int ycell)
{
	bin2D_t key_cell;
	key_cell.x = xcell;
	key_cell.y = ycell;
	return nb_container_exist(bins, &key_cell);
}

static inline uint32_t count_bins_in_block
                            (const nb_container_t *const restrict bins,
			     const int block[4])
{
	uint32_t N_bins;
	if (nb_container_get_length(bins) < 
	    (block[2] - block[0]) * (block[3] - block[1]))
		N_bins = count_bins_iterating_container(bins, block);
	else
		N_bins = count_bins_iterating_block(bins, block);
	return N_bins;
}

static uint32_t count_bins_iterating_container
                              (const nb_container_t *const restrict bins,
			       const int block[4])
{
	uint32_t N_bins = 0;
	uint32_t iter_size = nb_iterator_get_memsize();
	nb_iterator_t* iter = nb_soft_allocate_mem(iter_size);
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, bins);
	while (nb_iterator_has_more(iter)) {
		const bin2D_t* bin = nb_iterator_get_next(iter);
		if (bin_is_inside_block(bin->x, bin->y, block))
			N_bins += 1;
	}
	nb_iterator_finish(iter);
	nb_soft_free_mem(iter_size, iter);
	return N_bins;
}

static uint32_t count_bins_iterating_block
                              (const nb_container_t *const restrict bins,
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

static void get_bins_block(nb_container_t* bins_block,
			   const nb_container_t *const restrict bins,
			   int add_block[4],
			   int substract_block[4]
			   /* NULL if not required */)
{
	uint32_t bins_length = nb_container_get_length(bins);
	uint32_t N_blocks = (add_block[2] - add_block[0]) *
		(add_block[3] - add_block[1]);

	if (bins_length < N_blocks)
		iterate_cells_in_container(bins_block, bins, add_block,
					   substract_block);
	else
		iterate_cells_in_block(bins_block, bins, add_block,
				       substract_block);
}

static void iterate_cells_in_container(nb_container_t* bins_block,
				       const nb_container_t *bins,
				       int add_block[4],
				       int substract_block[4]
				       /* NULL if not required */)
{
	uint32_t iter_size = nb_iterator_get_memsize();
	nb_iterator_t* iter = nb_soft_allocate_mem(iter_size);
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, bins);
	while (nb_iterator_has_more(iter)) {
		const bin2D_t* bin = nb_iterator_get_next(iter);
		if (bin_is_outside_block(bin->x, bin->y, add_block))
			continue;
		if (NULL != substract_block)
			if (bin_is_inside_block(bin->x, bin->y,
						substract_block)) 
				continue;

		nb_container_insert(bins_block, bin->points);	
	}
	nb_iterator_finish(iter);
	nb_soft_free_mem(iter_size, iter);
}

static void iterate_cells_in_block(nb_container_t* bins_block,
				   const nb_container_t *bins,
				   int add_block[4],
				   int substract_block[4]
				   /* NULL if not required */)
{
	for (int i = add_block[0]; i <= add_block[2]; i++) {
		for (int j = add_block[1]; j <= add_block[3]; j++) {
			if (NULL != substract_block)
				if (bin_is_inside_block(i, j, substract_block)) 
					continue;
			bin2D_t *bin = get_bin(bins, i, j);
			if (NULL == bin)
				continue;

			nb_container_insert(bins_block, bin->points);
		}
	}
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

static void knn_get_layer_bins(const nb_bins2D_t *const restrict bins2D, 
			       int layer,
			       const nb_point2D_t *const restrict p,
			       nb_container_t *layer_bins)
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


	if (0 == layer) {
		bin2D_t* bin = get_bin(bins2D->bins, xbin, ybin);
		if (NULL != bin)
			/* Insert points */
			nb_container_insert(layer_bins, bin->points);
	} else if(0 < layer) {
		/* Get bottom horizontal bins */
		for (int i = xbin - layer; i <= xbin + layer; i++) {
			bin2D_t* bin = get_bin(bins2D->bins, i, ybin - layer);
			if (NULL != bin)
				nb_container_insert(layer_bins, bin->points);
		}
		/* Get top horizontal bins */
		for (int i = xbin - layer; i <= xbin + layer; i++) {
			bin2D_t* bin = get_bin(bins2D->bins, i, ybin + layer);
			if (NULL != bin)
				nb_container_insert(layer_bins, bin->points);
		}
		/* Get left vertical bins */
		for (int j = ybin - layer + 1; j < ybin + layer; j++) {
			bin2D_t* bin = get_bin(bins2D->bins, xbin - layer, j);
			if (NULL != bin)
				nb_container_insert(layer_bins, bin->points);
		}
		/* Get right vertical bins */
		for (int j = ybin - layer + 1; j < ybin + layer; j++) {
			bin2D_t* bin = get_bin(bins2D->bins, xbin + layer, j);
			if (NULL != bin)
				nb_container_insert(layer_bins, bin->points);
		}
	}
}

static inline void knn_scan_list(const nb_container_t *const points,
				 const nb_point2D_t *const restrict p_ref,
				 uint32_t k, nb_point2D_t** knn,
				 uint32_t* m_nearest,
				 double* knn_dist,
				 bool (*filter)(const nb_point2D_t *const p_ref,
						const nb_point2D_t *const p,
						const void *const data),
				 const void* const restrict filter_data)
{
	/* Iterate points in the cell */
	uint32_t iter_size = nb_iterator_get_memsize();
	nb_iterator_t* iter = nb_soft_allocate_mem(iter_size);
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, points);
	while (nb_iterator_has_more(iter)) {
		/* Get next point */
		nb_point2D_t *restrict p =
			(nb_point2D_t*) nb_iterator_get_next(iter);

		/* Filter if the pointer is the same */
		if (p_ref == p)
			continue;
    
		/* Filter point */
		if (!filter(p_ref, p, filter_data)) 
			continue;
    
		/* Calculate distance */
		double dist_p = nb_utils2D_get_dist(p_ref->x, p->x);

		for (uint32_t i = 0; i < m_nearest[0]; i++) {
			/* Minimize distance inserting in ascending order */
			if (dist_p < knn_dist[i]) {
				/* Swap value */
				double aux = dist_p;
				dist_p = knn_dist[i];
				knn_dist[i] = aux;
				nb_point2D_t* p_aux = p;
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
	nb_iterator_finish(iter);
	nb_soft_free_mem(iter_size, iter);
}

uint32_t nb_bins2D_get_knn(const nb_bins2D_t *const restrict bins2D,
			    const nb_point2D_t *const restrict p,
			    uint32_t k, nb_point2D_t* knn[], double knn_dist[])
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
	nb_container_type cnt_type = NB_QUEUE;
	uint32_t cnt_size = nb_container_get_memsize(cnt_type);
	nb_container_t *layer_cells = nb_soft_allocate_mem(cnt_size);

	uint32_t n_scanned = 0;
	uint32_t m_nearest = 0;
	uint32_t layer = 0;
	while (n_scanned < bins2D->length  && m_nearest < k) {
		/* 1. Get cells where the search should be performed. */
		nb_container_init(layer_cells, cnt_type);

		knn_get_layer_bins(bins2D, layer, p, layer_cells);

		/* 2. Compute nearest neighbours */
		while (nb_container_is_not_empty(layer_cells)) {
			nb_container_t *points =
				nb_container_delete_first(layer_cells);
			n_scanned += nb_container_get_length(points);
			knn_scan_list(points, p, k, knn, &m_nearest, knn_dist, 
				      bins2D->filter,
				      bins2D->filter_data);
		}
		/* Destroy list that had the layer cells */
		nb_container_finish(layer_cells);

		/* Increment layer */
		layer += 1;
	}
	nb_soft_free_mem(cnt_size, layer_cells);

	return m_nearest;
}

void nb_bins2D_set_filter(nb_bins2D_t *bins2D, 
			   bool (*filter)
			   (const nb_point2D_t *const p_ref,
			    const nb_point2D_t *const p,
			    const void *const data))
{
	bins2D->filter = filter;
}

void nb_bins2D_set_filter_data(nb_bins2D_t *bins2D, const void *data)
{
	bins2D->filter_data = data;
}

static bool half_side_have_points(const nb_bins2D_t* const bins2D,
				  const nb_point2D_t *const p1,
				  const nb_point2D_t *const p2)
{
	uint32_t iter_size = nb_iterator_get_memsize();
	nb_iterator_t* iter = nb_soft_allocate_mem(iter_size);
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, bins2D->bins);
	bool have_points = false;
	while (nb_iterator_has_more(iter)) {
		const bin2D_t* bin = nb_iterator_get_next(iter);
		if (nb_container_is_not_empty(bin->points)) {
			int8_t half_side = 
				bin_which_half_side(bin, bins2D->size_of_bins,
						    p1, p2);
			if (1 == half_side) {
				have_points = true;
				break;
			} else if (0 == half_side) {
				if (bin_have_points_in_half_side(bin, p1, p2)) {
					have_points = true;
					break;
				}
			}
		}
	}
	nb_iterator_finish(iter);
	nb_soft_free_mem(iter_size, iter);
	return have_points;
}

static int8_t bin_which_half_side(const bin2D_t* const bin,
				  double bin_size,
				  const nb_point2D_t *const p1,
				  const nb_point2D_t *const p2)
{
        double x1[2];
	x1[0] = bin->x * bin_size;
	x1[1] = bin->y * bin_size;

	double x2[2], x3[2], x4[2];
	x2[0] = x1[0] + bin_size;
	x2[1] = x1[1];
	x3[0] = x1[0] + bin_size;
	x3[1] = x1[1] + bin_size;
	x4[0] = x1[0];
	x4[1] = x1[1] + bin_size;

	bool h1 = nb_utils2D_is_in_half_side(p1->x, p2->x, x1);
	bool h2 = nb_utils2D_is_in_half_side(p1->x, p2->x, x2);
	bool h3 = nb_utils2D_is_in_half_side(p1->x, p2->x, x3);
	bool h4 = nb_utils2D_is_in_half_side(p1->x, p2->x, x4);

	int8_t half_side;
	if (h1 && h2 && h3 && h4)
		half_side = 1;
	else if (!h1 && !h2 && !h3 && !h4)
		half_side = -1;
	else
		half_side = 0;
	return half_side;
}

static bool bin_have_points_in_half_side(const bin2D_t* const bin,
					 const nb_point2D_t *const p1,
					 const nb_point2D_t *const p2)
{
	uint32_t iter_size = nb_iterator_get_memsize();
	nb_iterator_t* iter = nb_soft_allocate_mem(iter_size);
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, bin->points);
	bool have_points = false;
	while (nb_iterator_has_more(iter)) {
		const nb_point2D_t *p = nb_iterator_get_next(iter);
		if (nb_utils2D_is_in_half_side(p1->x, p2->x, p->x)) {
			have_points = true;
			break;
		}
	}
	nb_iterator_finish(iter);
	nb_soft_free_mem(iter_size, iter);
	return have_points;
}

static void delaunay_get_layer_bins
			(nb_container_t* layer_bins,
			 const nb_bins2D_t *const bins2D,
			 const nb_point2D_t *const p1,
			 const nb_point2D_t *const p2,
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
	get_bins_block(layer_bins, bins2D->bins,
		       layer_block, substract_block);
}

static void set_block_from_circle(int block[4],
				  const nb_bins2D_t *const bins2D,
				  const circle_t *const circle)
{
	block[0] = get_bin_coord(circle->c[0] - circle->r, bins2D->size_of_bins);
	block[1] = get_bin_coord(circle->c[1] - circle->r, bins2D->size_of_bins);
	block[2] = get_bin_coord(circle->c[0] + circle->r, bins2D->size_of_bins);
	block[3] = get_bin_coord(circle->c[1] + circle->r, bins2D->size_of_bins);

}

static void delaunay_get_points(const nb_container_t *const restrict points,
				const nb_point2D_t *const restrict p1,
				const nb_point2D_t *const restrict p2,
				const circle_t *const circle,
				nb_container_t* vertices,
				nb_container_t* outside_vtx)
{
	uint32_t iter_size = nb_iterator_get_memsize();
	nb_iterator_t* iter = nb_soft_allocate_mem(iter_size);
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, points);
	while (nb_iterator_has_more(iter)) {
		const nb_point2D_t* restrict p = nb_iterator_get_next(iter);
		delaunay_insert_if_is_candidate(vertices, outside_vtx,
						p1, p2, circle, p);
	}
	nb_iterator_finish(iter);
	nb_soft_free_mem(iter_size, iter);
}

static void delaunay_insert_if_is_candidate
				(nb_container_t* vertices,
				 nb_container_t* outside_vtx,
				 const nb_point2D_t *const restrict p1,
				 const nb_point2D_t *const restrict p2,
				 const circle_t *const circle,
				 const nb_point2D_t *const restrict p)
{
	if (p1 != p && p2 != p) {
		if (nb_utils2D_is_in_half_side(p1->x, p2->x, p->x)) {
			if (is_inside_circle(circle, p))
				nb_container_insert(vertices, p);
			else
				nb_container_insert(outside_vtx, p);
		}
	}
}

static inline bool is_inside_circle(const circle_t *const circle,
				    const nb_point2D_t *const p)
{
	double dist2 = nb_utils2D_get_dist2(circle->c, p->x);
	return (POW2(circle->r) > dist2);
}

void nb_bins2D_get_candidate_points_to_min_delaunay
                   (const nb_bins2D_t *const restrict bins2D,
		    const nb_point2D_t *const restrict p1, 
		    const nb_point2D_t *const restrict p2,
		    nb_container_t *vertices)
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
	nb_container_type cnt_type = NB_QUEUE;
	uint32_t cnt_size = nb_container_get_memsize(cnt_type);
	uint32_t memsize = 2 * cnt_size;
	char *memblock = nb_soft_allocate_mem(memsize);

	nb_container_t *outside_vtx = (void*) memblock;
	nb_container_t *layer_bins = (void*) (memblock + cnt_size);

	nb_container_init(outside_vtx, cnt_type);
	bool have_points = half_side_have_points(bins2D, p1, p2);
	uint32_t layer = 0;
	while (nb_container_is_empty(vertices) && have_points) {
		circle_t circle;
		set_circle_from_layer(&circle, p1, p2, bins2D, layer);

		delaunay_set_inside_points_of_prev_layer(vertices, outside_vtx,
							 &circle);

		/* 1. Get bins where the search should be performed. */
		nb_container_init(layer_bins, cnt_type);
		delaunay_get_layer_bins(layer_bins, bins2D, p1, p2,
					&circle, layer);

		/* 2. Get points from the bins of the layer */
		while (nb_container_is_not_empty(layer_bins)) {
			const nb_container_t *const restrict points =
				nb_container_delete_first(layer_bins);
			delaunay_get_points(points, p1, p2, &circle,
					    vertices, outside_vtx);
		}
		nb_container_finish(layer_bins);	
		layer += 1;

		if (layer > 200) { /**** TEMPORAL PATCH *************/
			fprintf(stderr, "-- NB_DEWALL WARNING: "	\
				"Bins layer search > 200, in "		\
				"Delaunay min dist\n");
			break;
		}/********************** TEMPORAL PATCH *************/
	}
	nb_container_finish(outside_vtx);

	nb_soft_free_mem(memsize, memblock);
}

static void set_circle_from_layer(circle_t *circle,
				  const nb_point2D_t *const p1,
				  const nb_point2D_t *const p2,
				  const nb_bins2D_t *const bins2D,
				  int layer)
{
	double n[2];
	n[0] = p1->x[1] - p2->x[1];
	n[1] = p2->x[0] - p1->x[0];

	double dist = sqrt(POW2(n[0]) + POW2(n[1]));
	n[0] /= dist;
	n[1] /= dist;

	double discret_dist = MAX(dist, bins2D->size_of_bins);
	circle->r = dist + layer * discret_dist;
	double a = sqrt(POW2(circle->r) - 0.25 * POW2(dist));
	circle->c[0] = 0.5 * (p1->x[0] + p2->x[0]) + a * n[0];
	circle->c[1] = 0.5 * (p1->x[1] + p2->x[1]) + a * n[1];
}

static void delaunay_set_inside_points_of_prev_layer
				(nb_container_t* vertices,
				 nb_container_t* outside_vtx,
				 const circle_t *const circle)
{
	nb_container_type type = nb_container_get_type(outside_vtx);
	uint32_t cnt_size = nb_container_get_memsize(type);
	nb_container_t *aux = nb_soft_allocate_mem(cnt_size);
	nb_container_init(aux, type);
	while (nb_container_is_not_empty(outside_vtx)) {
		nb_point2D_t *p = nb_container_delete_first(outside_vtx);
		if (is_inside_circle(circle, p))
			nb_container_insert(vertices, p);
		else
			nb_container_insert(aux, p);
	}
	nb_container_merge(outside_vtx, aux);
	nb_container_finish(aux);
	nb_soft_free_mem(cnt_size, aux);
}

void nb_bins2D_get_points_inside_circle(const nb_bins2D_t *const bins2D,
					double center[2],
					double radius,
					nb_container_t *points_inside)
{
	/* Get block */
	int block[4];
	block[0] = get_bin_coord(center[0] - radius, bins2D->size_of_bins);
	block[1] = get_bin_coord(center[1] - radius, bins2D->size_of_bins);
	block[2] = get_bin_coord(center[0] + radius, bins2D->size_of_bins);
	block[3] = get_bin_coord(center[1] + radius, bins2D->size_of_bins);

	nb_container_type cnt_type = NB_QUEUE;
	uint32_t cnt_size = nb_container_get_memsize(cnt_type);
	nb_container_t* bins_block = nb_soft_allocate_mem(cnt_size);
	nb_container_init(bins_block, cnt_type);

	get_bins_block(bins_block, bins2D->bins, block, NULL);

	get_points_inside_circle(bins_block, center, radius, points_inside);

	nb_container_finish(bins_block);
	nb_soft_free_mem(cnt_size, bins_block);
}

static void get_points_inside_circle(nb_container_t *bins_block,
				     double center[2], double radius,
				     nb_container_t *points_inside)
{
	uint32_t iter_size = nb_iterator_get_memsize();
	nb_iterator_t* iter = nb_soft_allocate_mem(iter_size);

	while (nb_container_is_not_empty(bins_block)) {
		nb_container_t* points =  nb_container_delete_first(bins_block);
		nb_iterator_init(iter);
		nb_iterator_set_container(iter, points);
		while (nb_iterator_has_more(iter)) {
			const nb_point2D_t *p = nb_iterator_get_next(iter);
			double dist2 = nb_utils2D_get_dist2(p->x, center);
			if (POW2(radius) - dist2 > NB_GEOMETRIC_TOL_POW2)
				nb_container_insert(points_inside, p);
		}
		nb_iterator_finish(iter);
	}

	nb_soft_free_mem(iter_size, iter);
}

bool nb_bins2D_are_points_inside_circle
                        (const nb_bins2D_t *const bins2D,
			 double center[2], double radius)
{
	/* Get block */
	int block[4];
	block[0]  = get_bin_coord(center[0] - radius, bins2D->size_of_bins);
	block[1]  = get_bin_coord(center[1] - radius, bins2D->size_of_bins);
	block[2]  = get_bin_coord(center[0] + radius, bins2D->size_of_bins);
	block[3]  = get_bin_coord(center[1] + radius, bins2D->size_of_bins);

	bool are_inside = false;
	for (int i = block[0]; i <= block[2]; i++) {
		for (int j = block[1]; j <= block[3]; j++) {
			bin2D_t *bin = get_bin(bins2D->bins, i, j);
			if (NULL == bin)
				continue;
			nb_container_t* points = bin->points;
			are_inside = are_points_inside_circle(points, center,
							      radius);
			if (are_inside)
				goto EXIT;
		}
	}
EXIT:
	return are_inside;
}

static bool are_points_inside_circle(const nb_container_t *points,
				     double center[2], double radius)
{
	uint32_t iter_size = nb_iterator_get_memsize();
	nb_iterator_t* iter = nb_soft_allocate_mem(iter_size);
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, points);
	bool are_inside = false;
	while (nb_iterator_has_more(iter)) {
		const nb_point2D_t *p = nb_iterator_get_next(iter);
		double dist2 = nb_utils2D_get_dist2(p->x, center);
		if (POW2(radius) - dist2 > NB_GEOMETRIC_TOL_POW2) {
			nb_iterator_finish(iter);
			are_inside = true;
			goto EXIT;
		}
	}
EXIT:
	nb_iterator_finish(iter);
	nb_soft_free_mem(iter_size, iter);
	return are_inside;
}

inline uint32_t nb_bins2D_get_N_bins(const nb_bins2D_t *const restrict bins2D)
{
	return nb_container_get_length(bins2D->bins);
}

inline bool nb_bins2D_is_empty(const nb_bins2D_t *const restrict bins2D)
{
	return nb_container_is_empty(bins2D->bins);
}

inline bool nb_bins2D_is_not_empty(const nb_bins2D_t *const restrict bins2D)
{
	return nb_container_is_not_empty(bins2D->bins);
}

uint32_t nb_bins2D_get_min_points_x_bin(const nb_bins2D_t *const bins2D)
{
	uint32_t iter_size = nb_iterator_get_memsize();
	nb_iterator_t* iter = nb_soft_allocate_mem(iter_size);
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, bins2D->bins);
	uint32_t min = bins2D->length;
	while (nb_iterator_has_more(iter)) {
		const bin2D_t *bin = nb_iterator_get_next(iter);
		if (min > nb_container_get_length(bin->points))
			min = nb_container_get_length(bin->points);
	}
	nb_iterator_finish(iter);
	nb_soft_free_mem(iter_size, iter);
	return min;
}

inline uint32_t nb_bins2D_get_length(const nb_bins2D_t *const restrict bins2D)
{
	return bins2D->length;
}

inline double nb_bins2D_get_size_of_bins
				(const nb_bins2D_t *const restrict bins2D)
{
	return bins2D->size_of_bins;
}
