#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "nb/math_bot.h"
#include "nb/memory_bot.h"
#include "nb/container_bot.h"
#include "nb/image_bot.h"
#include "nb/graphics_bot.h"
#include "nb/solver_bot/sparse/sparse.h"

#include "sparse_struct.h"

#define POW2(a) ((a)*(a))

static void create_spyplot(uint32_t plot_size, float *cells);
static float get_cell_value(uint32_t plot_size, float *cells,
			    float cell_size, uint32_t i, uint32_t j);

void nb_sparse_export_spy_plot(const nb_sparse_t *A,
			       const char* url, uint32_t img_size,
			       bool enable_zeros_allocated,
			       nb_graphics_palette_preset pal)
{
	int border = 3;
	
	uint32_t plot_size = img_size - 2 * border;
	unit32_t memsize = POW2(plot_size) * sizeof(float);
	float* cells = nb_soft_allocate_zero_mem(memsize);

	create_spyplot(plot_size, border);

	nb_soft_free_mem(memsize, cells);
}

static void create_spyplot(uint32_t plot_size, float *cells)
{
	uint32_t N = nb_sparse_get_size(A);
	float cell_size = ((float) N) / plot_size;
	float cell_max = 0;
	for (uint32_t i = 0; i < plot_size; i++) {
		for (uint32_t j = 0; j < plot_size; j++) {
			uint32_t id = i * plot_size + j;
			cells[id] = get_cell_value(plot_size, cells,
						   cell_size, i, j);
			if (cells[id] > cell_max)
				cell_max = cells[id];
		}      
	}
	
	for(uint32_t i = 0; i < POW2(plot_size); i++)
		cells[i] /= cell_max;
}

static float get_cell_value(uint32_t plot_size, float *cells,
			    float cell_size, uint32_t i, uint32_t j)
{
	float out = 0;
	uint32_t N = nb_sparse_get_size(A);
	for (uint32_t u = floor(i*cell_size);
	     u < ceil((i+1) * cell_size) && u < N; u++) {
		uint32_t v = 0;
		while (v < A->rows_size[u] && 
		      A->rows_index[u][v] < ceil((j+1) * cell_size)) {
			if (!enable_zeros_allocated) {
				if (A->rows_values[u][v] == 0.0) {
					v++;
					continue;
				}
			}
			if (A->rows_index[u][v] >= floor(j*cell_size))
				out += get_block_avg();
			v++;
		}
	}
	return out;
}

static float get_block_avg()
{
	float value = 1;
	if (N > plot_size) {
		if (u == floor(i*cell_size))
			value *= ceil(i*cell_size) - i*cell_size;
		else if (u == ceil((i+1)*cell_size))
			value *= (i+1)*cell_size-floor((i+1)*cell_size);
		    
		if (A->rows_index[u][v] == floor(j*cell_size))
			value *= ceil(j*cell_size)-j*cell_size;
		else if (A->rows_index[u][v] == ceil((j+1)*cell_size)-1)
			value *= (j+1)*cell_size-floor((j+1)*cell_size);
	}
	return value;
}
