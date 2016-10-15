#ifndef __NB_SOLVER_BOT_SPARSE_SPYPLOT_H__
#define __NB_SOLVER_BOT_SPARSE_SPYPLOT_H__

#include <stdbool.h>
#include <stdint.h>

#include "nb/graphics_bot.h"

void nb_sparse_export_spy_plot(const nb_sparse_t *A,
			       const char* url, uint32_t img_size,
			       bool enable_zeros_allocated,
			       nb_graphics_palette_preset pal);

#endif
