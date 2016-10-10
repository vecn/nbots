#ifndef __NB_GEOMETRIC_BOT_KNN_BINS2D_ITERATOR_H__
#define __NB_GEOMETRIC_BOT_KNN_BINS2D_ITERATOR_H__

#include <stdint.h>
#include <stdbool.h>
#include "nb/geometric_bot/point2D.h"

typedef struct nb_bins2D_iter_s nb_bins2D_iter_t;

uint32_t nb_bins2D_iter_get_memsize(void);
void nb_bins2D_iter_init(nb_bins2D_iter_t *iter);
void nb_bins2D_iter_finish(nb_bins2D_iter_t *iter);

nb_bins2D_iter_t* nb_bins2D_iter_create(void);
void nb_bins2D_iter_destroy(nb_bins2D_iter_t *iter);
void nb_bins2D_iter_set_bins(nb_bins2D_iter_t *iter,
			      const nb_bins2D_t *const bins2D);
bool nb_bins2D_iter_has_more(nb_bins2D_iter_t *iter);
const nb_point2D_t* nb_bins2D_iter_get_next(nb_bins2D_iter_t *iter);
void nb_bins2D_iter_restart(nb_bins2D_iter_t *iter);

#endif
