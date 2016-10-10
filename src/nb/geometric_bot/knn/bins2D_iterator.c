#include <stdlib.h>

#include "nb/container_bot.h"
#include "nb/geometric_bot/point2D.h"
#include "nb/geometric_bot/knn/bins2D.h"
#include "nb/geometric_bot/knn/bins2D_iterator.h"

#include "bins2D_structs.h"

struct nb_bins2D_iter_s {
	nb_container_t *bins;
	nb_iterator_t *bin_iter;
	nb_iterator_t *point_iter;
	const nb_point2D_t *current;
};

static void point_to_next(nb_bins2D_iter_t *iter);

uint32_t nb_bins2D_iter_get_memsize(void)
{
	return sizeof(nb_bins2D_iter_t);
}

void nb_bins2D_iter_init(nb_bins2D_iter_t *iter)
{
	memset(iter, 0, nb_bins2D_iter_get_memsize());
}

void nb_bins2D_iter_finish(nb_bins2D_iter_t *iter)
{
	if (NULL != iter->bin_iter)
		nb_iterator_destroy(iter->bin_iter);
	if (NULL != iter->point_iter)
		nb_iterator_destroy(iter->point_iter);
}

nb_bins2D_iter_t* nb_bins2D_iter_create(void)
{
	nb_bins2D_iter_t *iter = nb_allocate_mem(nb_bins2D_iter_get_memsize());
	nb_bins2D_iter_init(iter);
	return iter;
}

void nb_bins2D_iter_destroy(nb_bins2D_iter_t *iter)
{
	nb_bins2D_iter_finish(iter);
	nb_free_mem(iter);
}

inline void nb_bins2D_iter_set_bins(nb_bins2D_iter_t *iter,
				     const nb_bins2D_t *const bins2D)
{
	if (NULL != bins2D) {
		iter->bins = bins2D->bins;
		nb_bins2D_iter_restart(iter);
	}
}

inline bool nb_bins2D_iter_has_more(nb_bins2D_iter_t *iter)
{
	return (NULL != iter->current);
}

const nb_point2D_t* nb_bins2D_iter_get_next(nb_bins2D_iter_t *iter)
{
	const nb_point2D_t* point = iter->current;
	if (nb_iterator_has_more(iter->point_iter))
		iter->current = nb_iterator_get_next(iter->point_iter);
        else
		point_to_next(iter);
	return point;
}

static void point_to_next(nb_bins2D_iter_t *iter)
{
	iter->current = NULL;
	if (nb_iterator_has_more(iter->bin_iter)) {
		const bin2D_t* bin = nb_iterator_get_next(iter->bin_iter);
		nb_iterator_set_container(iter->point_iter, bin->points);
		if (nb_iterator_has_more(iter->point_iter))
			iter->current = nb_iterator_get_next(iter->point_iter);
	}
}

void nb_bins2D_iter_restart(nb_bins2D_iter_t *iter)
{
	if (NULL == iter->bin_iter)
		iter->bin_iter = nb_iterator_create();
	if (NULL == iter->point_iter)
		iter->point_iter = nb_iterator_create();

	nb_iterator_set_container(iter->bin_iter, iter->bins);
	point_to_next(iter);
}
