#ifndef __NB_SOLVER_BOT_SPARSE_SPYPLOT_H__
#define __NB_SOLVER_BOT_SPARSE_SPYPLOT_H__

#include <stdbool.h>
#include <stdint.h>

int vcn_sparse_spy_plot_as_png(const vcn_sparse_t *const A,
			       const char* url, uint32_t size,
			       bool enable_zeros_allocated,
			       bool enable_color);

#endif
