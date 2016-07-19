#ifndef __NB_GEOMETRIC_BOT_KNN_BINS2D_ITERATOR_H__
#define __NB_GEOMETRIC_BOT_KNN_BINS2D_ITERATOR_H__

#include <stdint.h>
#include <stdbool.h>
#include "nb/geometric_bot/point2D.h"

typedef struct vcn_bins2D_iter_s vcn_bins2D_iter_t;

uint32_t vcn_bins2D_iter_get_memsize(void);
void vcn_bins2D_iter_init(vcn_bins2D_iter_t *iter);
void vcn_bins2D_iter_finish(vcn_bins2D_iter_t *iter);

vcn_bins2D_iter_t* vcn_bins2D_iter_create(void);
void vcn_bins2D_iter_destroy(vcn_bins2D_iter_t *iter);
void vcn_bins2D_iter_set_bins(vcn_bins2D_iter_t *iter,
			      const vcn_bins2D_t *const bins2D);
bool vcn_bins2D_iter_has_more(vcn_bins2D_iter_t *iter);
const vcn_point2D_t* vcn_bins2D_iter_get_next(vcn_bins2D_iter_t *iter);
void vcn_bins2D_iter_restart(vcn_bins2D_iter_t *iter);

#endif
