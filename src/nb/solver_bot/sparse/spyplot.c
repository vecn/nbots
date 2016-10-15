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
#include "nb/solver_bot/sparse/spyplot.h"

#include "sparse_struct.h"

#define POW2(a) ((a)*(a))


static void init_image(nb_image_t *img, uint32_t size, uint8_t *pix);
static void create_spyplot(const nb_sparse_t *A,
			   uint32_t plot_size, float *cells,
			   bool enable_zeros_allocated);
static float get_cell_value(const nb_sparse_t *A,
			    uint32_t plot_size, float *cells,
			    float cell_size, uint32_t i, uint32_t j,
			    bool enable_zeros_allocated);
static float get_block_avg(const nb_sparse_t *A,
			   uint32_t plot_size, float *cells,
			   float cell_size, uint32_t i, uint32_t j,
			   uint32_t u, uint32_t v);
static void set_image(nb_image_t *img, int border, uint32_t plot_size,
		      float *cells, nb_palette_preset pal);

void nb_sparse_export_spy_plot(const nb_sparse_t *A,
			       const char* url, uint32_t img_size,
			       bool enable_zeros_allocated,
			       nb_palette_preset pal)
{
	int border = 3;
	
	uint32_t plot_size = img_size - 2 * border;
	
	uint32_t img_memsize = nb_image_get_memsize();
	uint32_t pix_memsize = POW2(plot_size) * 4 * sizeof(uint8_t);
	uint32_t memsize = POW2(plot_size) * sizeof(float) +
		img_memsize + pix_memsize;
	char *memblock = nb_soft_allocate_mem(memsize);
	nb_image_t *img = (void*) memblock;
	uint8_t *pix = (void*) (memblock  + img_memsize);
	float* cells =  (void*) (memblock + img_memsize + pix_memsize);
	
	init_image(img, img_size, pix);

	create_spyplot(A, plot_size, cells, enable_zeros_allocated);

	set_image(img, border, plot_size, cells, pal);

	nb_image_write(img, url);

	nb_soft_free_mem(memsize, memblock);
}

static void init_image(nb_image_t *img, uint32_t size, uint8_t *pix)
{
	nb_image_init(img);
	img->comp_x_pixel = 4;
	img->width = size;
	img->height = size;
	img->pixels = pix;
}

static void create_spyplot(const nb_sparse_t *A,
			   uint32_t plot_size, float *cells,
			   bool enable_zeros_allocated)
{
	memset(cells, 0, POW2(plot_size) * sizeof(*cells));
	uint32_t N = nb_sparse_get_size(A);
	float cell_size = ((float) N) / plot_size;
	float cell_max = 0;
	for (uint32_t i = 0; i < plot_size; i++) {
		for (uint32_t j = 0; j < plot_size; j++) {
			uint32_t id = i * plot_size + j;
			cells[id] = get_cell_value(A, plot_size, cells,
						   cell_size, i, j,
						   enable_zeros_allocated);
			if (cells[id] > cell_max)
				cell_max = cells[id];
		}      
	}
	
	for(uint32_t i = 0; i < POW2(plot_size); i++)
		cells[i] /= cell_max;
}

static float get_cell_value(const nb_sparse_t *A,
			    uint32_t plot_size, float *cells,
			    float cell_size, uint32_t i, uint32_t j,
			    bool enable_zeros_allocated)
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
				out += get_block_avg(A, plot_size, cells,
						     cell_size,
						     i, j, u, v);
			v++;
		}
	}
	return out;
}

static float get_block_avg(const nb_sparse_t *A,
			   uint32_t plot_size, float *cells,
			   float cell_size, uint32_t i, uint32_t j,
			   uint32_t u, uint32_t v)
{
	float value = 1;
	if (A->N > plot_size) {
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

static void set_image(nb_image_t *img, int border, uint32_t plot_size,
		      float *cells, nb_palette_preset pal)
{
	nb_palette_t *palette =
		nb_palette_create_preset(pal);
	
	for (uint32_t i = 0; i < img->width; i++) {
		for (uint32_t j = 0; j < img->height; j++) {
			uint8_t pix[4] = {0, 0, 0, 255};
			if (i >= border && i <  border + plot_size &&
			    j >= border && j <  border + plot_size) {
				uint32_t r = i - border;
				uint32_t c = j - border;
				uint32_t id = r *  plot_size + c;
				nb_palette_get_rgba(palette,
							     cells[id],
							     pix);
			}
			nb_image_set_pixel(img, i, j, pix);
		}
	}

	nb_palette_destroy(palette);
}
